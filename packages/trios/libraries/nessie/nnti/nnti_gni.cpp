/**
 * nnti_gni.c
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
#include <ifaddrs.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netdb.h>

#include <sched.h>
#include <alps/libalpslli.h>
#include <gni_pub.h>

#include <map>
#include <deque>

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

#include "nnti_gni.h"
#include "nnti_utils.h"



///* if undefined, the ACK message is NOT sent to the RDMA target when
// * the RDMA op is complete.  this creates one-side semantics for RDMA
// * ops.  in this mode, the target has no idea when the RDMA op is
// * complete and what data was addressed.
// */
//#undef USE_RDMA_TARGET_ACK
///* if defined, the RDMA initiator will send an ACK message to the RDMA
// * target when the RDMA op is complete.  the target process must wait
// * on the target buffer in order to get the ACK.  this creates two-side
// * semantics for RDMA ops.   in this mode, when the wait returns the
// * the RDMA op is complete and status indicates what data was addressed.
// */
//#define USE_RDMA_TARGET_ACK



//#undef USE_ANONYMOUS_MEMORY_REGISTRATION
//#define USE_ANONYMOUS_MEMORY_REGISTRATION

//#undef USE_RDMA_EVENTS
//#define USE_RDMA_EVENTS

// NSSI
//#define USE_RDMA_TARGET_ACK
//#define USE_ANONYMOUS_MEMORY_REGISTRATION
//#define USE_RDMA_EVENTS
//#undef USE_RDMA_FENCE

// ORNL
//#undef USE_RDMA_TARGET_ACK
//#define USE_ANONYMOUS_MEMORY_REGISTRATION
//#undef USE_RDMA_EVENTS
//#undef USE_RDMA_FENCE


//#define USE_FMA
//#define USE_BTE
//#define USE_CROSSOVER
//#define USE_MIXED

//#define FMA_BTE_CROSSOVER 4096

///* mode 1 - client uses GNI params from ALPS to create a Commumnication Domain and attach to it */
//#define USE_ALPS_PTAG
///* mode 2 - client uses a mix of GNI params from ALPS and the server (cookie/pTag) to create a Commumnication Domain and attach to it */
////#undef USE_ALPS_PTAG



#define NIC_ADDR_BITS    22
#define NIC_ADDR_SHIFT   (32-NIC_ADDR_BITS)
#define NIC_ADDR_MASK    0x3FFFFF
#define CPU_NUM_BITS     7
#define CPU_NUM_SHIFT    (NIC_ADDR_SHIFT-CPU_NUM_BITS)
#define CPU_NUM_MASK     0x7F
#define THREAD_NUM_BITS  3
#define THREAD_NUM_SHIFT (CPU_NUM_SHIFT-THREAD_NUM_BITS)
#define THREAD_NUM_MASK  0x7

#define GNI_INSTID(nic_addr, cpu_num, thr_num) (((nic_addr&NIC_ADDR_MASK)<<NIC_ADDR_SHIFT)|((cpu_num&CPU_NUM_MASK)<<CPU_NUM_SHIFT)|(thr_num&THREAD_NUM_MASK))
#define GNI_NIC_ADDRESS(inst_id)               ((inst_id>>NIC_ADDR_SHIFT)&NIC_ADDR_MASK)
#define GNI_CPU_NUMBER(inst_id)                ((inst_id>>CPU_NUM_SHIFT)&CPU_NUM_MASK)
#define GNI_THREAD_NUMBER(inst_id)             (inst_id&THREAD_NUM_MASK)



enum nnti_gni_pi_ordering {
    PI_ORDERING_DEFAULT,
    PI_ORDERING_STRICT,
    PI_ORDERING_RELAXED
};
enum nnti_gni_rdma_mode {
    RDMA_FMA,
    RDMA_BTE,
    RDMA_MIXED,
    RDMA_CROSSOVER
};
typedef struct {

    bool use_alps_ptag;

    bool use_wr_pool;

    bool use_rdma_target_ack;
    bool use_rdma_events;
    bool use_rdma_fence;

    enum nnti_gni_pi_ordering pi_ordering; /* DEFAULT, STRICT, RELAXED */
    enum nnti_gni_rdma_mode   rdma_mode;   /* FMA, BTE, MIXED, CROSSOVER */

    uint32_t fma_bte_crossover_size;

} nnti_gni_config;


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
} gni_connection_state;

typedef enum {
    SERVER_CONNECTION,
    CLIENT_CONNECTION
} gni_connection_type;

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
} gni_buffer_type;


#define GNI_OP_PUT_INITIATOR  1
#define GNI_OP_GET_INITIATOR  2
#define GNI_OP_PUT_TARGET     3
#define GNI_OP_GET_TARGET     4
#define GNI_OP_SEND_REQUEST   5
#define GNI_OP_SEND_BUFFER    6
#define GNI_OP_NEW_REQUEST    7
#define GNI_OP_RECEIVE        8

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
} gni_op_state_t;

#define SRQ_WQ_DEPTH 2048
#define CQ_WQ_DEPTH 128
#define ACK_PER_CONN 64

typedef struct {
    gni_cq_handle_t cq_hdl;
    gni_ep_handle_t ep_hdl;
} conn_ep;


typedef struct {
    uint64_t ack_received;
    uint64_t inst_id;
    uint64_t byte_len;
    uint64_t byte_offset;
    uint64_t src_offset;
    uint64_t dest_offset;
    uint16_t op;
} nnti_gni_work_completion;

/**
 * attrs to send to the client
 */
typedef struct {
    uint64_t          req_index_addr; /* address of the request index var on the server */
    gni_mem_handle_t  req_index_mem_hdl;

    uint64_t         req_buffer_addr; /* address of the request buffer */
    uint64_t         req_size;        /* the maximum size of a request */
    uint64_t         req_count;       /* the number of requests that will fit in the queue */
    gni_mem_handle_t req_mem_hdl;

    uint64_t         wc_buffer_addr;  /* address of the work completion array */
    gni_mem_handle_t wc_mem_hdl;
} nnti_gni_recv_queue_attrs;

/**
 * attrs to send to the server
 */
typedef struct {
    uint64_t         unblock_buffer_addr; /* address of the unblock buffer */
    gni_mem_handle_t unblock_mem_hdl;     /* mem hdl to send to the server */
} nnti_gni_recv_queue_flow_control_attrs;

typedef struct {
    nnti_gni_recv_queue_flow_control_attrs client;
    nnti_gni_recv_queue_attrs              server;
} nnti_gni_queue_remote_attrs;


typedef struct {
    uint64_t         req_index;      /* after AMO Fetch-Add, this buffer contains the current index of the request queue */
    uint64_t         req_index_addr; /* address of the request index var on the client. */
    gni_cq_handle_t  req_index_mem_cq_hdl;
    gni_mem_handle_t req_index_mem_hdl;
    gni_cq_handle_t  req_index_cq_hdl;
    gni_ep_handle_t  req_index_ep_hdl;

    gni_cq_handle_t  req_cq_hdl;
    gni_ep_handle_t  req_ep_hdl;

    uint64_t         unblock_buffer;      /* small buffer to receive unblock messages */
    uint64_t         unblock_buffer_addr; /* address of the unblock buffer */
    gni_cq_handle_t  unblock_mem_cq_hdl;  /* CQ to wait for unblock events */
    gni_mem_handle_t unblock_mem_hdl;     /* mem hdl to send to the server.  the server uses this mem hdl as the
                                           * remote hdl when posting (CQ write) the unblock message. */

    uint64_t last_offset;
} nnti_gni_client_queue;

typedef struct {
    NNTI_peer_t       peer;
    char             *peer_name;
    NNTI_ip_addr      peer_addr;
    NNTI_tcp_port     peer_port;
    uint32_t          peer_cookie;
    uint32_t          peer_ptag;
    NNTI_instance_id  peer_instance;

    alpsAppGni_t      peer_alps_info;

    nnti_gni_client_queue       queue_local_attrs;
    nnti_gni_queue_remote_attrs queue_remote_attrs;

    gni_cdm_handle_t     cdm_hdl;  /* a client creates this comm domain with params from the server */
    gni_nic_handle_t     nic_hdl;
    gni_ep_handle_t      ep_hdl;

    gni_connection_state state;

    gni_connection_type connection_type;
} gni_connection;

typedef struct {
    uint8_t                  is_initiator;

    const NNTI_buffer_t     *reg_buf;

    gni_post_descriptor_t    post_desc;
    gni_post_descriptor_t   *post_desc_ptr;

    nnti_gni_work_completion wc;
    gni_mem_handle_t         wc_mem_hdl;
    uint64_t                 wc_dest_addr;
    gni_mem_handle_t         wc_dest_mem_hdl;
    gni_post_descriptor_t    wc_post_desc;
    int8_t                   wc_registered;

    uint8_t                  last_op;
    uint8_t                  is_op_complete;
    gni_op_state_t           op_state;

    NNTI_instance_id         peer_instance;
} gni_work_request;

typedef std::deque<gni_work_request *>           wr_queue_t;
typedef std::deque<gni_work_request *>::iterator wr_queue_iter_t;

typedef struct {
    gni_buffer_type  type;
    gni_mem_handle_t mem_hdl;
    wr_queue_t       wr_queue;
    nthread_lock_t   wr_queue_lock;
    uint32_t         ref_count;
} gni_memory_handle;

typedef struct {
    NNTI_buffer_t            *reg_buf;

    uint64_t                  last_index_before_reset;
    uint64_t                  req_index;       /* index of the next available slot in the request queue */
    uint64_t                  req_index_addr;  /* address of the request index var on the server */
    gni_cq_handle_t           req_index_mem_cq_hdl;
    gni_mem_handle_t          req_index_mem_hdl;

    char                     *req_buffer;      /* pointer to the head of the request buffer */
    uint64_t                  req_size;        /* the maximum size of a request */
    uint64_t                  req_count;       /* the number of requests that will fit in the queue */
    uint64_t                  req_buffer_size; /* the size of the request buffer in bytes (req_size*req_count) */

    nnti_gni_work_completion *wc_buffer;       /* pointer to the head of the work completion array */
    uint64_t                  wc_buffer_size;  /* the size of the work completion buffer in bytes (sizeof(nnti_gni_work_completion)*req_count) */
    gni_cq_handle_t           wc_mem_cq_hdl;
    gni_mem_handle_t          wc_mem_hdl;

    gni_cq_handle_t           unblock_cq_hdl;
    gni_ep_handle_t           unblock_ep_hdl;

    uint64_t                  req_processed_reset_limit;
    uint64_t                  req_processed;
    uint64_t                  total_req_processed;
} gni_request_queue_handle;

typedef struct {
    uint16_t         delivery_mode;

    gni_cdm_handle_t cdm_hdl;
    gni_nic_handle_t nic_hdl;
    NNTI_instance_id instance;

    gni_cq_handle_t  ep_cq_hdl;
    gni_cq_handle_t  mem_cq_hdl;
    gni_cq_handle_t  cq_list[2];

    uint64_t         apid;
    alpsAppGni_t     alps_info;

    int              listen_sock;
    char             listen_name[NNTI_HOSTNAME_LEN];
    uint32_t         listen_addr;  /* in NBO */
    uint16_t         listen_port;  /* in NBO */

    gni_request_queue_handle req_queue;
} gni_transport_global;




static nthread_lock_t nnti_gni_lock;


static NNTI_result_t register_memory(
        gni_memory_handle *hdl,
        void *buf,
        uint64_t len);
static NNTI_result_t unregister_memory(
        gni_memory_handle *hdl);
static NNTI_result_t register_wc(
        gni_work_request *wr);
static NNTI_result_t unregister_wc(
        gni_work_request *wr);
static void reset_op_state(
        const NNTI_buffer_t *reg_buf);
static gni_cq_handle_t get_cq(
        const NNTI_buffer_t *reg_buf);
static const NNTI_buffer_t *decode_event_buffer(
        const NNTI_buffer_t *wait_buf,
        gni_cq_entry_t      *ev_data);
static int process_event(
        const NNTI_buffer_t      *reg_buf,
        gni_cq_handle_t           cq_hdl,
        gni_cq_entry_t           *ev_data);
static NNTI_result_t post_recv_work_request(
        NNTI_buffer_t *reg_buf);
static NNTI_result_t repost_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        gni_work_request *wr);
static gni_work_request *first_incomplete_wr(
        gni_memory_handle *gni_mem_hdl);
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
        int                   nnti_rc,
        gni_cq_entry_t       *ev_data,
        NNTI_status_t        *status);
static void create_peer(NNTI_peer_t *peer,
        char *name,
        NNTI_ip_addr addr,
        NNTI_tcp_port port,
        uint32_t ptag,
        uint32_t cookie,
        NNTI_instance_id instance);
static void copy_peer(
        NNTI_peer_t *src,
        NNTI_peer_t *dest);
static int init_server_listen_socket(void);
static int check_listen_socket_for_new_connections(void);
static uint32_t get_cpunum(void);
static void get_alps_info(alpsAppGni_t *alps_info);
static int tcp_read(int sock, void *incoming, size_t len);
static int tcp_write(int sock, const void *outgoing, size_t len);
static int tcp_exchange(int sock, int is_server, void *incoming, void *outgoing, size_t len);
static void transition_connection_to_ready(
        int sock,
        gni_connection *conn);
static NNTI_result_t get_ipaddr(
        char *ipaddr,
        int maxlen);
static NNTI_result_t init_connection(
        gni_connection **conn,
        const int sock,
        const int is_server);
static void close_connection(gni_connection *c);
static void print_wc(
        const nnti_gni_work_completion *wc);
static void print_cq_event(
        const gni_cq_entry_t *event);
static void print_post_desc(
        const gni_post_descriptor_t *post_desc_ptr);
static int need_mem_cq(gni_memory_handle *gni_mem_hdl);
static int need_wc_mem_cq(gni_memory_handle *gni_mem_hdl);
static void print_failed_cq(const NNTI_buffer_t *reg_buf);
static void print_gni_conn(gni_connection *c);
static NNTI_result_t insert_conn_peer(const NNTI_peer_t *peer, gni_connection *conn);
static NNTI_result_t insert_conn_instance(const NNTI_instance_id instance, gni_connection *conn);
static gni_connection *get_conn_peer(const NNTI_peer_t *peer);
static gni_connection *get_conn_instance(const NNTI_instance_id instance);
static NNTI_peer_t *get_peer_by_url(const char *url);
static gni_connection *del_conn_peer(const NNTI_peer_t *peer);
static gni_connection *del_conn_instance(const NNTI_instance_id instance);
static void print_peer_map(void);
static void print_instance_map(void);
static void close_all_conn(void);
//static void print_put_buf(void *buf, uint32_t size);
static void print_raw_buf(void *buf, uint32_t size);

static NNTI_result_t insert_buf_bufhash(NNTI_buffer_t *buf);
static NNTI_buffer_t *get_buf_bufhash(const uint32_t bufhash);
static NNTI_buffer_t *del_buf_bufhash(NNTI_buffer_t *buf);
static void print_bufhash_map(void);

static NNTI_result_t insert_wr_wrhash(gni_work_request *);
static gni_work_request *get_wr_wrhash(const uint32_t bufhash);
static gni_work_request *del_wr_wrhash(gni_work_request *);
static void print_wrhash_map(void);

static NNTI_result_t wr_pool_init(uint32_t pool_size);
static gni_work_request *wr_pool_target_pop(void);
static gni_work_request *wr_pool_initiator_pop(void);
static void wr_pool_target_push(gni_work_request *wr);
static void wr_pool_initiator_push(gni_work_request *wr);
static NNTI_result_t wr_pool_fini(void);

static uint16_t get_dlvr_mode_from_env();
static void set_dlvr_mode(
        gni_post_descriptor_t *pd);
static void set_rdma_mode(
        gni_post_descriptor_t *pd);
static void set_post_desc(
        gni_post_descriptor_t *pd,
        gni_buffer_type target_buf_type,
        uint32_t buf_length);
static void set_wc_post_desc(
        gni_post_descriptor_t *pd,
        uint32_t buf_length);

static int server_req_queue_init(
        gni_request_queue_handle *q,
        char                     *buffer,
        uint64_t                  req_size,
        uint64_t                  req_count);
static int server_req_queue_destroy(
        gni_request_queue_handle *q);

static int client_req_queue_init(
        gni_connection *c);
static int client_req_queue_destroy(
        gni_connection *c);

static int send_unblock(
        gni_request_queue_handle *local_req_queue_attrs);
static int reset_req_index(
        gni_request_queue_handle  *req_queue_attrs);

static void send_rdma_wc (
        gni_work_request    *wr,
        const NNTI_buffer_t *local_buf,
        const NNTI_buffer_t *remote_buf);

static int fetch_add_buffer_offset(
        nnti_gni_client_queue       *local_req_queue_attrs,
        nnti_gni_recv_queue_attrs *remote_req_queue_attrs,
        uint64_t                     addend,
        uint64_t                    *prev_offset);
static int send_req(
        const NNTI_peer_t           *peer_hdl,
        nnti_gni_client_queue       *local_req_queue_attrs,
        nnti_gni_recv_queue_attrs *remote_req_queue_attrs,
        uint64_t                     offset,
        const NNTI_buffer_t         *reg_buf,
        gni_work_request            *wr);
static int send_req_wc(
        const NNTI_peer_t           *peer_hdl,
        nnti_gni_client_queue       *local_req_queue_attrs,
        nnti_gni_recv_queue_attrs *remote_req_queue_attrs,
        uint64_t                     offset,
        const NNTI_buffer_t         *reg_buf,
        gni_work_request            *wr);
static int request_send(
        const NNTI_peer_t           *peer_hdl,
        nnti_gni_client_queue       *client_q,
        nnti_gni_recv_queue_attrs *server_q,
        const NNTI_buffer_t         *reg_buf,
        int                          req_num);
static int send_buffer(
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *src_hdl,
        const NNTI_buffer_t *dest_hdl,
        gni_work_request    *wr);
static int send_buffer_wc(
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *src_hdl,
        const NNTI_buffer_t *dest_hdl,
        gni_work_request    *wr);
static int buffer_send(
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *src_hdl,
        const NNTI_buffer_t *dest_hdl);

static void config_init(
        nnti_gni_config *c);
static void config_get_from_env(
        nnti_gni_config *c);


static bool gni_initialized=false;

static gni_transport_global transport_global_data;
static const int MIN_TIMEOUT = 1;  /* in milliseconds.  must be >0 for GNI_CqVectorWaitEvent(). */

static log_level nnti_cq_debug_level;
static log_level nnti_event_debug_level;
static log_level nnti_ee_debug_level;



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
static std::map<addrport_key, gni_connection *> connections_by_peer;
typedef std::map<addrport_key, gni_connection *>::iterator conn_by_peer_iter_t;
typedef std::pair<addrport_key, gni_connection *> conn_by_peer_t;
static nthread_lock_t nnti_conn_peer_lock;

static std::map<NNTI_instance_id, gni_connection *> connections_by_instance;
typedef std::map<NNTI_instance_id, gni_connection *>::iterator conn_by_inst_iter_t;
typedef std::pair<NNTI_instance_id, gni_connection *> conn_by_inst_t;
static nthread_lock_t nnti_conn_instance_lock;

static std::map<uint32_t, NNTI_buffer_t *> buffers_by_bufhash;
typedef std::map<uint32_t, NNTI_buffer_t *>::iterator buf_by_bufhash_iter_t;
typedef std::pair<uint32_t, NNTI_buffer_t *> buf_by_bufhash_t;
static nthread_lock_t nnti_buf_bufhash_lock;

static std::map<uint32_t, gni_work_request *> wr_by_wrhash;
typedef std::map<uint32_t, gni_work_request *>::iterator wr_by_wrhash_iter_t;
typedef std::pair<uint32_t, gni_work_request *> wr_by_wrhash_t;
static nthread_lock_t nnti_wr_wrhash_lock;

typedef std::deque<gni_work_request *>           wr_pool_t;
typedef std::deque<gni_work_request *>::iterator wr_pool_iter_t;
static nthread_lock_t nnti_wr_pool_lock;

static wr_pool_t target_wr_pool;
static wr_pool_t initiator_wr_pool;


static nnti_gni_config config;


/* ---------------- Wrappers to protect HPCToolkit issues ---------- */
static
gni_return_t GNI_CqWaitEvent_wrapper(
        gni_cq_handle_t cq_hdl,
        uint64_t timeout_per_call,
        gni_cq_entry_t *ev_data)
{
    gni_return_t  rc;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc = GNI_CqWaitEvent(cq_hdl, timeout_per_call, ev_data);
    if (sampling) SAMPLING_START();

    return rc;
}

static
gni_return_t GNI_CdmAttach_wrapper(
        gni_cdm_handle_t cdm_hndl,
        uint32_t device_id,
        uint32_t *local_addr,
        gni_nic_handle_t *nic_hndl)
{
    gni_return_t rc;
    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=GNI_CdmAttach (cdm_hndl, device_id, local_addr, nic_hndl);
    if (sampling) SAMPLING_START();
    return rc;
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
NNTI_result_t NNTI_gni_init (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl)
{
    int rc=NNTI_OK;

    trios_declare_timer(call_time);

    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
//    char memdesc[NNTI_URL_LEN];
    char *sep, *endptr;

    char hostname[NNTI_HOSTNAME_LEN];
    NNTI_ip_addr  addr;
    NNTI_tcp_port port;

    uint32_t nic_addr  =0;
    uint32_t cpu_num   =0;
    uint32_t thread_num=0;
    uint32_t gni_cpu_id=0;


    assert(trans_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    log_debug(nnti_debug_level, "my_url=%s", my_url);
    log_debug(nnti_debug_level, "initialized=%d, FALSE==%d", (int)gni_initialized, (int)FALSE);

    if (!gni_initialized) {

        memset(&transport_global_data, 0, sizeof(gni_transport_global));

        nnti_cq_debug_level=nnti_debug_level;
        nnti_event_debug_level=nnti_debug_level;
        nnti_ee_debug_level=nnti_debug_level;

        nthread_lock_init(&nnti_gni_lock);

        // initialize the mutexes for the connection maps
        nthread_lock_init(&nnti_conn_peer_lock);
        nthread_lock_init(&nnti_conn_instance_lock);
        nthread_lock_init(&nnti_wr_wrhash_lock);
        nthread_lock_init(&nnti_buf_bufhash_lock);

        nthread_lock_init(&nnti_wr_pool_lock);

        config_init(&config);
        config_get_from_env(&config);

        log_debug(nnti_debug_level, "my_url=%s", my_url);

        hostname[0]='\0';
        if (my_url != NULL) {
            if ((rc=nnti_url_get_transport(my_url, transport, NNTI_URL_LEN)) != NNTI_OK) {
                goto cleanup;
            }
            if (0!=strcmp(transport, "gni")) {
                rc=NNTI_EINVAL;
                goto cleanup;
            }

            if ((rc=nnti_url_get_address(my_url, address, NNTI_URL_LEN)) != NNTI_OK) {
                goto cleanup;
            }

//            if ((rc=nnti_url_get_memdesc(my_url, memdesc, NNTI_URL_LEN)) != NNTI_OK) {
//                return(rc);
//            }

            sep=strchr(address, ':');
            if (sep == address) {
                /* no hostname given; try gethostname */
                gethostname(hostname, NNTI_HOSTNAME_LEN);
            } else {
                strncpy(hostname, address, sep-address);
                hostname[sep-address]='\0';
            }
            sep++;
            port=strtol(sep, &endptr, 0);
            if (endptr == sep) {
                /* no port given; use -1 */
                port=-1;
            }
        } else {
            rc=get_ipaddr(hostname, NNTI_HOSTNAME_LEN);
            if (rc != NNTI_OK) {
                log_error(nnti_debug_level, "could not find IP address to listen on");
                goto cleanup;
            }
//            gethostname(hostname, NNTI_HOSTNAME_LEN);
            port=-1;
        }
        strcpy(transport_global_data.listen_name, hostname);


        transport_global_data.delivery_mode = get_dlvr_mode_from_env();


        log_debug(nnti_debug_level, "initializing Gemini");

//        /* register trace groups (let someone else enable) */
//        trace_counter_gid = trace_get_gid(TRACE_RPC_COUNTER_GNAME);
//        trace_interval_gid = trace_get_gid(TRACE_RPC_INTERVAL_GNAME);


        trios_start_timer(call_time);
        get_alps_info(&transport_global_data.alps_info);
        trios_stop_timer("get_alps_info", call_time);

        trios_start_timer(call_time);
        rc=GNI_CdmGetNicAddress (transport_global_data.alps_info.device_id, &nic_addr, &gni_cpu_id);
        trios_stop_timer("CdmGetNicAddress", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CdmGetNicAddress() failed: %d", rc);
            if (rc==GNI_RC_NO_MATCH)
                rc=NNTI_EEXIST;
            else
                rc=NNTI_EINVAL;

            goto cleanup;
        }

        trios_start_timer(call_time);
        cpu_num = get_cpunum();
        trios_stop_timer("get_cpunum", call_time);

        transport_global_data.instance=GNI_INSTID(nic_addr, cpu_num, thread_num);
        log_debug(nnti_debug_level, "nic_addr(%llu), cpu_num(%llu), thread_num(%llu), inst_id(%llu), "
                "derived.nic_addr(%llu), derived.cpu_num(%llu), derived.thread_num(%llu)",
                (uint64_t)nic_addr, (uint64_t)cpu_num, (uint64_t)thread_num,
                (uint64_t)transport_global_data.instance,
                (uint64_t)GNI_NIC_ADDRESS(transport_global_data.instance),
                (uint64_t)GNI_CPU_NUMBER(transport_global_data.instance),
                (uint64_t)GNI_THREAD_NUMBER(transport_global_data.instance));

        log_debug(nnti_debug_level, "global_nic_hdl - host(%s) device_id(%llu) local_addr(%lld) cookie(%llu) ptag(%llu) "
                    "apid(%llu) inst_id(%llu) gni_nic_addr(%llu) gni_cpu_id(%llu) linux_cpu_num(%llu) omp_thread_num(%llu)",
                    transport_global_data.listen_name,
                    (unsigned long long)transport_global_data.alps_info.device_id,
                    (long long)transport_global_data.alps_info.local_addr,
                    (unsigned long long)transport_global_data.alps_info.cookie,
                    (unsigned long long)transport_global_data.alps_info.ptag,
                    (unsigned long long)transport_global_data.apid,
                    (unsigned long long)transport_global_data.instance,
                    (unsigned long long)nic_addr,
                    (unsigned long long)gni_cpu_id,
                    (unsigned long long)cpu_num,
                    (unsigned long long)thread_num);

        trios_start_timer(call_time);
        rc=GNI_CdmCreate(transport_global_data.instance,
                transport_global_data.alps_info.ptag,
                transport_global_data.alps_info.cookie,
                GNI_CDM_MODE_ERR_NO_KILL,
                &transport_global_data.cdm_hdl);
        trios_stop_timer("CdmCreate", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CdmCreate() failed: %d", rc);
            rc=NNTI_EINVAL;
            goto cleanup;
        }

        trios_start_timer(call_time);

        rc=GNI_CdmAttach_wrapper (transport_global_data.cdm_hdl,
                transport_global_data.alps_info.device_id,
                (uint32_t*)&transport_global_data.alps_info.local_addr, /* ALPS and GNI disagree about the type of local_addr.  cast here. */
                &transport_global_data.nic_hdl);

        trios_stop_timer("CdmAttach", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CdmAttach() failed: %d", rc);
            if (rc==GNI_RC_PERMISSION_ERROR)
                rc=NNTI_EPERM;
            else
                rc=NNTI_EINVAL;

            goto cleanup;
        }

        rc=GNI_CqCreate (transport_global_data.nic_hdl, 1000, 0, GNI_CQ_BLOCKING, NULL, NULL, &transport_global_data.ep_cq_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CqCreate(transport_global_data.ep_cq_hdl) failed: %d", rc);
            goto cleanup;
        }
        rc=GNI_CqCreate (transport_global_data.nic_hdl, 1000, 0, GNI_CQ_BLOCKING, NULL, NULL, &transport_global_data.mem_cq_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CqCreate(transport_global_data.mem_cq_hdl) failed: %d", rc);
            goto cleanup;
        }
        transport_global_data.cq_list[0]=transport_global_data.ep_cq_hdl;
        transport_global_data.cq_list[1]=transport_global_data.mem_cq_hdl;

        if (config.use_wr_pool) {
            rc=wr_pool_init(100);
            if (rc!=NNTI_OK) {
                log_error(nnti_debug_level, "wr_pool_init(): %d", rc);
                goto cleanup;
            }
        }

        trios_start_timer(call_time);
        init_server_listen_socket();
        trios_stop_timer("init_server_listen_socket", call_time);
//        trios_start_timer(call_time);
//        start_connection_listener_thread();
//        trios_stop_timer("start_connection_listener_thread", call_time);

        if (logging_info(nnti_debug_level)) {
            fprintf(logger_get_file(), "Gemini Initialized: host(%s) port(%u)\n",
                    transport_global_data.listen_name,
                    ntohs(transport_global_data.listen_port));
        }

        trios_start_timer(call_time);
        create_peer(
                &trans_hdl->me,
                transport_global_data.listen_name,
                transport_global_data.listen_addr,
                transport_global_data.listen_port,
                transport_global_data.alps_info.ptag,
                transport_global_data.alps_info.cookie,
                transport_global_data.instance);
        trios_stop_timer("create_peer", call_time);

        gni_initialized = true;
    }

cleanup:
    log_debug(nnti_ee_debug_level, "exit");

    return((NNTI_result_t)rc);
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
 *                - ex. "ptl://nid:pid/", "gni://ip_addr:port", "luc://endpoint_id/"
 *    - memory_descriptor - (optional) transport-specific representation of RMA params
 */
NNTI_result_t NNTI_gni_get_url (
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen)
{
    NNTI_result_t rc=NNTI_OK;

    assert(trans_hdl);
    assert(url);
    assert(maxlen>0);

    log_debug(nnti_ee_debug_level, "enter");

    strncpy(url, trans_hdl->me.url, maxlen);
    url[maxlen-1]='\0';

    log_debug(nnti_ee_debug_level, "exit");

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
NNTI_result_t NNTI_gni_connect (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl)
{
    int rc=NNTI_OK;

    trios_declare_timer(call_time);

    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
    char params[NNTI_URL_LEN];
    char *sep;

    NNTI_peer_t *existing_peer=NULL;

    char          hostname[NNTI_HOSTNAME_LEN];
    char          port_str[NNTI_HOSTNAME_LEN];
    NNTI_tcp_port port;

    char     *cookie_str;
    char     *ptag_str;

    int s;
    struct hostent *host_entry;
    struct sockaddr_in skin;
    socklen_t skin_size=sizeof(struct sockaddr_in);

    gni_connection *conn=NULL;

    NNTI_peer_t *key;

    long start_time;
    long elapsed_time = 0;
    long timeout_per_call;


    assert(trans_hdl);
    assert(peer_hdl);

    log_debug(nnti_ee_debug_level, "enter (url=%s)", url);

    existing_peer=get_peer_by_url(url);
    if (existing_peer!=NULL) {
        *peer_hdl=*existing_peer;
        return(NNTI_OK);
    }

    conn = (gni_connection *)calloc(1, sizeof(gni_connection));
    log_debug(nnti_debug_level, "calloc returned conn=%p.", conn);
    if (conn == NULL) {
        log_error(nnti_debug_level, "calloc returned NULL.  out of memory?: %s", strerror(errno));
        rc=NNTI_ENOMEM;
        goto cleanup;
    }

    if (url != NULL) {
        if ((rc=nnti_url_get_transport(url, transport, NNTI_URL_LEN)) != NNTI_OK) {
            goto cleanup;
        }
        if (0!=strcmp(transport, "gni")) {
            /* the peer described by 'url' is not an Gemini peer */
            rc=NNTI_EINVAL;
            goto cleanup;
        }

        if ((rc=nnti_url_get_address(url, address, NNTI_URL_LEN)) != NNTI_OK) {
            /* the peer described by 'url' is not an Gemini peer */
            rc=NNTI_EINVAL;
            goto cleanup;
        }

        if ((rc=nnti_url_get_params(url, params, NNTI_URL_LEN)) != NNTI_OK) {
            /* the peer described by 'url' is not an Gemini peer */
            rc=NNTI_EINVAL;
            goto cleanup;
        }

        sep=strchr(address, ':');
        strncpy(hostname, address, sep-address);
        hostname[sep-address]='\0';
        strcpy(port_str, sep+1);
        port=strtol(port_str, NULL, 0);

        log_debug(nnti_ee_debug_level, "params=%s", params);

        ptag_str=strstr(params, "ptag=");
        sep=strchr(ptag_str, '&');
        *sep='\0';
        log_debug(nnti_ee_debug_level, "ptag_str=%s", ptag_str+5);
        conn->peer_ptag=strtol(ptag_str+5, NULL, 10);
        *sep='&';

        cookie_str=strstr(params, "cookie=");
        log_debug(nnti_ee_debug_level, "cookie_str=%s", cookie_str+7);
        conn->peer_cookie=strtoull(cookie_str+7, NULL, 10);

        log_debug(nnti_ee_debug_level, "url=%s", url);

    } else {
        /*  */
        rc=NNTI_EINVAL;
        goto cleanup;
    }

    if (config.use_alps_ptag) {
        conn->cdm_hdl = transport_global_data.cdm_hdl;
        conn->nic_hdl = transport_global_data.nic_hdl;
    } else {
        log_debug(nnti_ee_debug_level, "conn_nic_hdl - host(%s) device_id(%llu) local_addr(%lld) cookie(%llu) ptag(%llu) "
                "apid(%llu) inst_id(%llu) gni_nic_addr(%llu) linux_cpu_num(%llu) omp_thread_num(%llu)",
                transport_global_data.listen_name,
                (unsigned long long)transport_global_data.alps_info.device_id,
                (long long)transport_global_data.alps_info.local_addr,
                (unsigned long long)conn->peer_cookie,
                (unsigned long long)conn->peer_ptag,
                (unsigned long long)transport_global_data.apid,
                (unsigned long long)transport_global_data.instance,
                (uint64_t)GNI_NIC_ADDRESS(transport_global_data.instance),
                (uint64_t)GNI_CPU_NUMBER(transport_global_data.instance),
                (uint64_t)GNI_THREAD_NUMBER(transport_global_data.instance));

        trios_start_timer(call_time);
        rc=GNI_CdmCreate(transport_global_data.instance,
                conn->peer_ptag,
                conn->peer_cookie,
                GNI_CDM_MODE_ERR_ALL_KILL,
                &conn->cdm_hdl);
        trios_stop_timer("CdmCreate", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CdmCreate() failed: %d", rc);
            rc=NNTI_EINVAL;
            goto cleanup;
        }

        trios_start_timer(call_time);
        rc=GNI_CdmAttach(conn->cdm_hdl,
                transport_global_data.alps_info.device_id,
                (uint32_t*)&transport_global_data.alps_info.local_addr, /* ALPS and GNI disagree about the type of local_addr.  cast here. */
                &conn->nic_hdl);
        trios_stop_timer("CdmAttach", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CdmAttach() failed: %d", rc);
            if (rc==GNI_RC_PERMISSION_ERROR)
                rc=NNTI_EPERM;
            else
                rc=NNTI_EINVAL;
            goto cleanup;
        }

        rc=GNI_CqCreate (conn->nic_hdl, 1, 0, GNI_CQ_BLOCKING, NULL, NULL, &transport_global_data.ep_cq_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CqCreate(transport_global_data.ep_cq_hdl) failed: %d", rc);
            goto cleanup;
        }
        rc=GNI_CqCreate (conn->nic_hdl, 1, 0, GNI_CQ_BLOCKING, NULL, NULL, &transport_global_data.mem_cq_hdl);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "CqCreate(transport_global_data.mem_cq_hdl) failed: %d", rc);
            goto cleanup;
        }
        transport_global_data.cq_list[0]=transport_global_data.ep_cq_hdl;
        transport_global_data.cq_list[1]=transport_global_data.mem_cq_hdl;
    }

    host_entry = gethostbyname(hostname);
    if (!host_entry) {
        log_warn(nnti_debug_level, "failed to resolve server name (%s): %s", hostname, strerror(errno));
        rc=NNTI_ENOENT;
        goto cleanup;
    }
    memset(&skin, 0, sizeof(skin));
    skin.sin_family = host_entry->h_addrtype;
    memcpy(&skin.sin_addr, host_entry->h_addr_list[0], (size_t) host_entry->h_length);
    skin.sin_port = htons(port);

    elapsed_time=0;
    timeout_per_call = MIN_TIMEOUT;

    s = socket(AF_INET, SOCK_STREAM, 0);
    if (s < 0) {
        log_warn(nnti_debug_level, "failed to create tcp socket: errno=%d (%s)", errno, strerror(errno));
        rc=NNTI_EIO;
        goto cleanup;
    }
    trios_start_timer(call_time);
    start_time=trios_get_time_ms();
    while((timeout==-1) || (elapsed_time < timeout)) {
        log_debug(nnti_debug_level, "calling connect");
        if (connect(s, (struct sockaddr *)&skin, skin_size) == 0) {
            log_debug(nnti_debug_level, "connected");
            break;
        }
        elapsed_time=trios_get_time_ms()-start_time;
        log_warn(nnti_debug_level, "failed to connect to server (%s:%u): errno=%d (%s)", hostname, port, errno, strerror(errno));
        if ((timeout>0) && (elapsed_time >= timeout)) {
            rc=NNTI_EIO;
            goto cleanup;
        }
        nnti_sleep(timeout_per_call);
    }
    trios_stop_timer("socket connect", call_time);

    trios_start_timer(call_time);
    rc = init_connection(&conn, s, 0);
    trios_stop_timer("gni init connection", call_time);
    if (conn==NULL) {
        rc=NNTI_EIO;
        goto cleanup;
    }

    create_peer(
            peer_hdl,
            conn->peer_name,
            conn->peer_addr,
            conn->peer_port,
            conn->peer_ptag,
            conn->peer_cookie,
            conn->peer_instance);

    conn->peer=*peer_hdl;

    insert_conn_peer(peer_hdl, conn);
    insert_conn_instance(conn->peer_instance, conn);

    transition_connection_to_ready(s, conn);

    if (close(s) < 0) {
        log_warn(nnti_debug_level, "failed to close tcp socket: errno=%d (%s)", errno, strerror(errno));
        rc=NNTI_EIO;
        goto cleanup;
    }

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "peer_hdl",
                "end of NNTI_gni_connect", peer_hdl);
    }

cleanup:
    if (rc != NNTI_OK) {
        if (conn!=NULL) free(conn);
    }
    log_debug(nnti_ee_debug_level, "exit");

    return((NNTI_result_t)rc);
}


/**
 * @brief Terminate communication with this peer.
 *
 * Perform any transport specific actions necessary to end communication with
 * this peer.
 */
NNTI_result_t NNTI_gni_disconnect (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    assert(trans_hdl);
    assert(peer_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    gni_connection *conn=get_conn_peer(peer_hdl);
    close_connection(conn);
    del_conn_peer(peer_hdl);

    log_debug(nnti_ee_debug_level, "exit");

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
NNTI_result_t NNTI_gni_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        const NNTI_peer_t      *peer,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;
    trios_declare_timer(call_time);

    NNTI_buffer_t     *old_buf=NULL;
    gni_memory_handle *gni_mem_hdl=NULL;

    assert(trans_hdl);
    assert(buffer);
    assert(element_size>0);
    assert(num_elements>0);
    assert(ops>0);
    assert(reg_buf);

    log_debug(nnti_ee_debug_level, "enter");

    old_buf=get_buf_bufhash(hash6432shift((uint64_t)buffer));
    if (old_buf==NULL) {
        gni_mem_hdl=new gni_memory_handle();
        assert(gni_mem_hdl);
        gni_mem_hdl->ref_count=1;
        nthread_lock_init(&gni_mem_hdl->wr_queue_lock);
    } else {
        gni_mem_hdl=(gni_memory_handle*)old_buf->transport_private;
        gni_mem_hdl->ref_count++;
    }

    log_debug(nnti_ee_debug_level, "gni_mem_hdl->ref_count==%lu", gni_mem_hdl->ref_count);

    reg_buf->transport_id      = trans_hdl->id;
    reg_buf->buffer_owner      = trans_hdl->me;
    reg_buf->ops               = ops;
    reg_buf->payload_size      = element_size;
    reg_buf->payload           = (uint64_t)buffer;
    reg_buf->transport_private = (uint64_t)gni_mem_hdl;

    log_debug(nnti_debug_level, "rpc_buffer->payload_size=%ld",
            reg_buf->payload_size);

    if (gni_mem_hdl->ref_count==1) {
        if (ops == NNTI_RECV_QUEUE) {
            gni_request_queue_handle *q_hdl=&transport_global_data.req_queue;

            gni_mem_hdl->type   =REQUEST_BUFFER;
            //        gni_mem_hdl->last_op=GNI_OP_NEW_REQUEST;

            memset(q_hdl, 0, sizeof(gni_request_queue_handle));

            q_hdl->reg_buf=reg_buf;

            server_req_queue_init(
                    q_hdl,
                    buffer,
                    element_size,
                    num_elements);

            reg_buf->payload_size=q_hdl->req_size;

        } else if (ops == NNTI_RECV_DST) {
            gni_mem_hdl->type    =RECEIVE_BUFFER;
        } else if (ops == NNTI_SEND_SRC) {
            gni_mem_hdl->type=SEND_BUFFER;
        } else if (ops == NNTI_GET_DST) {
            gni_mem_hdl->type    =GET_DST_BUFFER;
        } else if (ops == NNTI_GET_SRC) {
            gni_mem_hdl->type    =GET_SRC_BUFFER;
        } else if (ops == NNTI_PUT_SRC) {
            gni_mem_hdl->type    =PUT_SRC_BUFFER;
        } else if (ops == NNTI_PUT_DST) {
            gni_mem_hdl->type    =PUT_DST_BUFFER;
        } else if (ops == (NNTI_GET_SRC|NNTI_PUT_DST)) {
            gni_mem_hdl->type    =RDMA_TARGET_BUFFER;
        } else {
            gni_mem_hdl->type=UNKNOWN_BUFFER;
        }

        if (ops != NNTI_RECV_QUEUE) {
            rc=register_memory(gni_mem_hdl, buffer, element_size);
        }
    }

    if ((ops == NNTI_RECV_QUEUE) || (ops == NNTI_RECV_DST)) {
        post_recv_work_request(reg_buf);
    }

    if (config.use_rdma_target_ack) {
        if ((gni_mem_hdl->type == RDMA_TARGET_BUFFER) ||
                (gni_mem_hdl->type == GET_SRC_BUFFER) ||
                (gni_mem_hdl->type == PUT_DST_BUFFER)) {
            post_recv_work_request(reg_buf);
        }
    }

    if (rc==NNTI_OK) {
        reg_buf->buffer_addr.transport_id                            = NNTI_TRANSPORT_GEMINI;
        reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword1 = gni_mem_hdl->mem_hdl.qword1;
        reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword2 = gni_mem_hdl->mem_hdl.qword2;

        if (gni_mem_hdl->type==REQUEST_BUFFER) {
            reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.size = transport_global_data.req_queue.req_size;
            reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf  = (uint64_t)transport_global_data.req_queue.req_buffer;
            reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.type = NNTI_GNI_REQUEST_BUFFER;
        } else {
            reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.size = reg_buf->payload_size;
            reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf  = (uint64_t)reg_buf->payload;
            reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.type = NNTI_GNI_SEND_SRC;
        }
    }


    if (gni_mem_hdl->ref_count==1) {
        insert_buf_bufhash(reg_buf);
        log_debug(nnti_debug_level, "gni_mem_hdl->type==%llu",
                (uint64_t)gni_mem_hdl->type);
        log_debug(nnti_debug_level, "reg_buf.buf.hash==%llu",
                (uint64_t)hash6432shift(reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf));
    }

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_gni_register_memory", reg_buf);
    }

    log_debug(nnti_ee_debug_level, "exit");
    return(rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t NNTI_gni_unregister_memory (
        NNTI_buffer_t    *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;
    gni_memory_handle *gni_mem_hdl=NULL;

    assert(reg_buf);

    log_debug(nnti_ee_debug_level, "enter");

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "start of NNTI_gni_unregister_memory", reg_buf);
    }

    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
    assert(gni_mem_hdl);
    gni_mem_hdl->ref_count--;

    log_debug(nnti_ee_debug_level, "gni_mem_hdl->ref_count==%lu", gni_mem_hdl->ref_count);

    if (gni_mem_hdl->ref_count==0) {
        log_debug(nnti_ee_debug_level, "gni_mem_hdl->ref_count is 0.  release all resources.");

        if (gni_mem_hdl->type==REQUEST_BUFFER) {
            server_req_queue_destroy(
                    &transport_global_data.req_queue);

        } else {
            unregister_memory(gni_mem_hdl);
        }

        del_buf_bufhash(reg_buf);

        nthread_lock(&gni_mem_hdl->wr_queue_lock);
        while (!gni_mem_hdl->wr_queue.empty()) {
            gni_work_request *wr=gni_mem_hdl->wr_queue.front();
            log_debug(nnti_debug_level, "removing pending wr=%p", wr);
            gni_mem_hdl->wr_queue.pop_front();
            del_wr_wrhash(wr);
            if (config.use_wr_pool) {
                if (wr->is_initiator==TRUE) {
                    wr_pool_initiator_push(wr);
                } else {
                    wr_pool_target_push(wr);
                }
            } else {
                unregister_wc(wr);
                free(wr);
            }
        }
        nthread_unlock(&gni_mem_hdl->wr_queue_lock);

        nthread_lock_fini(&gni_mem_hdl->wr_queue_lock);

        if (gni_mem_hdl) delete gni_mem_hdl;

        reg_buf->transport_id      = NNTI_TRANSPORT_NULL;
        GNI_SET_MATCH_ANY(&reg_buf->buffer_owner);
        reg_buf->ops               = (NNTI_buf_ops_t)0;
        //    GNI_SET_MATCH_ANY(&reg_buf->peer);
        reg_buf->payload_size      = 0;
        reg_buf->payload           = 0;
        reg_buf->transport_private = 0;
    }

    log_debug(nnti_ee_debug_level, "exit");

    return(rc);
}


/**
 * @brief Send a message to a peer.
 *
 * Send a message (<tt>msg_hdl</tt>) to a peer (<tt>peer_hdl</tt>).  It is expected that the
 * message is small, but the exact maximum size is transport dependent.
 */
NNTI_result_t NNTI_gni_send (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    trios_declare_timer(call_time);

    gni_memory_handle *gni_mem_hdl=NULL;

    log_debug(nnti_ee_debug_level, "enter");

    assert(peer_hdl);
    assert(msg_hdl);

    gni_mem_hdl=(gni_memory_handle *)msg_hdl->transport_private;
    assert(gni_mem_hdl);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "msg_hdl",
                "NNTI_send", msg_hdl);
    }
    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "dest_hdl",
                "NNTI_send", dest_hdl);
    }

    if ((dest_hdl == NULL) || (dest_hdl->ops == NNTI_RECV_QUEUE)) {
        gni_connection *conn=get_conn_peer(peer_hdl);
        assert(conn);

        trios_start_timer(call_time);
        request_send(peer_hdl, &conn->queue_local_attrs, &conn->queue_remote_attrs.server, msg_hdl, 0);
        trios_stop_timer("send to request queue", call_time);

    } else if (dest_hdl->ops == NNTI_RECV_DST) {
        gni_connection *conn=get_conn_peer(peer_hdl);
        assert(conn);

        trios_start_timer(call_time);
        buffer_send(peer_hdl, msg_hdl, dest_hdl);
        trios_stop_timer("send to receive dest", call_time);
    }

    log_debug(nnti_debug_level, "sending to (%s, instance=%llu)", peer_hdl->url, (uint64_t)peer_hdl->peer.NNTI_remote_process_t_u.gni.inst_id);

    log_debug(nnti_ee_debug_level, "exit");

    return(rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * Put the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_gni_put (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset)
{
    int rc=NNTI_OK;
    trios_declare_timer(call_time);

    gni_ep_handle_t        ep_hdl;

    gni_memory_handle *gni_mem_hdl=NULL;
    gni_work_request  *wr=NULL;

    bool use_fma=true;

    log_debug(nnti_ee_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);


    gni_mem_hdl=(gni_memory_handle *)src_buffer_hdl->transport_private;
    assert(gni_mem_hdl);

    if (config.use_wr_pool) {
        wr=wr_pool_initiator_pop();
    } else {
        wr=(gni_work_request *)calloc(1, sizeof(gni_work_request));
    }
    assert(wr);

    wr->reg_buf = src_buffer_hdl;

    memset(&wr->post_desc, 0, sizeof(gni_post_descriptor_t));
    set_post_desc(&wr->post_desc, PUT_DST_BUFFER, src_length);

    wr->post_desc.local_addr            =src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf+src_offset;
    wr->post_desc.local_mem_hndl.qword1 =src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
    wr->post_desc.local_mem_hndl.qword2 =src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword2;
    wr->post_desc.remote_addr           =dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf+dest_offset;
    wr->post_desc.remote_mem_hndl.qword1=dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
    wr->post_desc.remote_mem_hndl.qword2=dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword2;
    wr->post_desc.length                =src_length;


    wr->op_state=RDMA_WRITE_INIT;

    if (!config.use_wr_pool) {
        wr->wc_registered  =FALSE;
    }
    wr->wc.op         =GNI_OP_PUT_TARGET;
    wr->wc.byte_len   =src_length;
    wr->wc.src_offset =src_offset;
    wr->wc.dest_offset=dest_offset;
    wr->wc.inst_id    =transport_global_data.instance;

    wr->wc_dest_addr          =dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.wc_addr;
    wr->wc_dest_mem_hdl.qword1=dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword1;
    wr->wc_dest_mem_hdl.qword2=dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword2;

    gni_connection *conn=get_conn_peer(&dest_buffer_hdl->buffer_owner);
    assert(conn);
    ep_hdl=conn->ep_hdl;

    nthread_lock(&gni_mem_hdl->wr_queue_lock);

    nthread_lock(&nnti_gni_lock);

    GNI_EpSetEventData(
            ep_hdl,
            hash6432shift((uint64_t)wr),
            hash6432shift((uint64_t)dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf));

    if (config.rdma_mode==RDMA_CROSSOVER) {
        if (wr->post_desc.length > config.fma_bte_crossover_size) {
            use_fma=false;
        } else {
            use_fma=true;
        }
    } else if ((config.rdma_mode==RDMA_BTE) || (config.rdma_mode==RDMA_MIXED)) {
        use_fma=false;
    } else if (config.rdma_mode==RDMA_FMA) {
        use_fma=true;
    }

    if (use_fma) {
        log_debug(nnti_event_debug_level, "calling PostFma(fma put ; ep_hdl(%llu) transport_global_data.cq_hdl(%llu) local_mem_hdl(%llu, %llu) remote_mem_hdl(%llu, %llu))",
                ep_hdl, transport_global_data.ep_cq_hdl,
                wr->post_desc.local_mem_hndl.qword1, wr->post_desc.local_mem_hndl.qword2,
                wr->post_desc.remote_mem_hndl.qword1, wr->post_desc.remote_mem_hndl.qword2);
        trios_start_timer(call_time);
        rc=GNI_PostFma(ep_hdl, &wr->post_desc);
        trios_stop_timer("PostFma put", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "failed to post FMA (rc=%d): %s", rc, strerror(errno));
            rc=NNTI_EIO;
        }
    } else {
        log_debug(nnti_event_debug_level, "calling PostRdma(rdma put ; ep_hdl(%llu) transport_global_data.cq_hdl(%llu) local_mem_hdl(%llu, %llu) remote_mem_hdl(%llu, %llu))",
                ep_hdl, transport_global_data.ep_cq_hdl,
                wr->post_desc.local_mem_hndl.qword1, wr->post_desc.local_mem_hndl.qword2,
                wr->post_desc.remote_mem_hndl.qword1, wr->post_desc.remote_mem_hndl.qword2);
        trios_start_timer(call_time);
        rc=GNI_PostRdma(ep_hdl, &wr->post_desc);
        trios_stop_timer("PostRdma put", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "failed to post RDMA (rc=%d): %s", rc, strerror(errno));
            rc=NNTI_EIO;
        }
    }

    wr->last_op  =GNI_OP_PUT_INITIATOR;
    wr->peer_instance=dest_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.gni.inst_id;
    log_debug(nnti_event_debug_level, "wr->peer_instance==%lu", wr->peer_instance);

    if (config.use_rdma_target_ack) {
        send_rdma_wc(wr, src_buffer_hdl, dest_buffer_hdl);
    }

    nthread_unlock(&nnti_gni_lock);

    gni_mem_hdl->wr_queue.push_back(wr);
    insert_wr_wrhash(wr);

    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

    log_debug(nnti_ee_debug_level, "exit");

    return((NNTI_result_t)rc);
}


/**
 * @brief Transfer data from a peer.
 *
 * Get the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_gni_get (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset)
{
    int rc=NNTI_OK;
    trios_declare_timer(call_time);

    gni_ep_handle_t        ep_hdl;

    gni_memory_handle *gni_mem_hdl=NULL;
    gni_work_request  *wr=NULL;

    bool use_fma=true;

    log_debug(nnti_ee_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);


    gni_mem_hdl=(gni_memory_handle *)dest_buffer_hdl->transport_private;
    assert(gni_mem_hdl);

    if (config.use_wr_pool) {
        wr=wr_pool_initiator_pop();
    } else {
        wr=(gni_work_request *)calloc(1, sizeof(gni_work_request));
    }
    assert(wr);

    wr->reg_buf = dest_buffer_hdl;

    memset(&wr->post_desc, 0, sizeof(gni_post_descriptor_t));
    set_post_desc(&wr->post_desc, GET_SRC_BUFFER, src_length);

    gni_connection *conn=get_conn_peer(&src_buffer_hdl->buffer_owner);
    assert(conn);
    ep_hdl=conn->ep_hdl;

    nthread_lock(&gni_mem_hdl->wr_queue_lock);

    nthread_lock(&nnti_gni_lock);

    GNI_EpSetEventData(
            ep_hdl,
            hash6432shift((uint64_t)wr),
            hash6432shift((uint64_t)src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf));

    wr->post_desc.local_addr            =dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf+dest_offset;
    wr->post_desc.local_mem_hndl.qword1 =dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
    wr->post_desc.local_mem_hndl.qword2 =dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword2;
    wr->post_desc.remote_addr           =src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf+src_offset;
    wr->post_desc.remote_mem_hndl.qword1=src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
    wr->post_desc.remote_mem_hndl.qword2=src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword2;
    wr->post_desc.length                =src_length;

    wr->op_state=RDMA_READ_INIT;

    if (!config.use_wr_pool) {
        wr->wc_registered  =FALSE;
    }
    wr->wc.op         =GNI_OP_GET_TARGET;
    wr->wc.byte_len   =src_length;
    wr->wc.src_offset =src_offset;
    wr->wc.dest_offset=dest_offset;
    wr->wc.inst_id    =transport_global_data.instance;

    wr->wc_dest_addr          =src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.wc_addr;
    wr->wc_dest_mem_hdl.qword1=src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword1;
    wr->wc_dest_mem_hdl.qword2=src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword2;

    if (config.rdma_mode==RDMA_CROSSOVER) {
        if (wr->post_desc.length > config.fma_bte_crossover_size) {
            use_fma=false;
        } else {
            use_fma=true;
        }
    } else if ((config.rdma_mode==RDMA_BTE) || (config.rdma_mode==RDMA_MIXED)) {
        use_fma=false;
    } else if (config.rdma_mode==RDMA_FMA) {
        use_fma=true;
    }

    if (use_fma) {
        log_debug(nnti_event_debug_level, "calling PostFma(fma get ; ep_hdl(%llu) transport_global_data.cq_hdl(%llu) local_mem_hdl(%llu, %llu) remote_mem_hdl(%llu, %llu))",
                ep_hdl, transport_global_data.ep_cq_hdl,
                wr->post_desc.local_mem_hndl.qword1, wr->post_desc.local_mem_hndl.qword2,
                wr->post_desc.remote_mem_hndl.qword1, wr->post_desc.remote_mem_hndl.qword2);
        trios_start_timer(call_time);
        rc=GNI_PostFma(ep_hdl, &wr->post_desc);
        trios_stop_timer("PostFma get", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "failed to post FMA (rc=%d): %s", rc, strerror(errno));
            rc=NNTI_EIO;
        }
    } else {
        log_debug(nnti_event_debug_level, "calling PostRdma(rdma get ; ep_hdl(%llu) transport_global_data.cq_hdl(%llu) local_mem_hdl(%llu, %llu) remote_mem_hdl(%llu, %llu))",
                ep_hdl, transport_global_data.ep_cq_hdl,
                wr->post_desc.local_mem_hndl.qword1, wr->post_desc.local_mem_hndl.qword2,
                wr->post_desc.remote_mem_hndl.qword1, wr->post_desc.remote_mem_hndl.qword2);
        trios_start_timer(call_time);
        rc=GNI_PostRdma(ep_hdl, &wr->post_desc);
        trios_stop_timer("PostRdma get", call_time);
        if (rc!=GNI_RC_SUCCESS) {
            log_error(nnti_debug_level, "failed to post RDMA (rc=%d): %s", rc, strerror(errno));
            rc=NNTI_EIO;
        }
    }

    wr->last_op  =GNI_OP_GET_INITIATOR;
    wr->peer_instance=src_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.gni.inst_id;
    log_debug(nnti_event_debug_level, "wr->peer_instance==%lu", wr->peer_instance);

    if (config.use_rdma_target_ack) {
        send_rdma_wc(wr, dest_buffer_hdl, src_buffer_hdl);
    }

    nthread_unlock(&nnti_gni_lock);

    gni_mem_hdl->wr_queue.push_back(wr);
    insert_wr_wrhash(wr);

    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

    log_debug(nnti_ee_debug_level, "exit");

    return((NNTI_result_t)rc);
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
NNTI_result_t NNTI_gni_wait (
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t        *status)
{
    int nnti_rc=NNTI_OK;
    gni_memory_handle        *gni_mem_hdl=NULL;
    gni_work_request         *wr=NULL;
    gni_request_queue_handle *q_hdl=NULL;
    gni_connection           *conn=NULL;
    const NNTI_buffer_t      *wait_buf=NULL;

    gni_ep_handle_t        ep_hdl;

    gni_cq_handle_t cq_hdl=0;
    gni_cq_entry_t  ev_data;

//    nnti_gni_work_completion wc;

    gni_return_t rc=GNI_RC_SUCCESS;
    long elapsed_time = 0;
    long timeout_per_call;

    long entry_time=trios_get_time_ms();

    trios_declare_timer(call_time);

    log_debug(nnti_ee_debug_level, "enter");

    assert(reg_buf);
    assert(status);

    q_hdl      =&transport_global_data.req_queue;
    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
    assert(q_hdl);
    assert(gni_mem_hdl);

    check_listen_socket_for_new_connections();

    if (!config.use_rdma_target_ack) {
        if ((remote_op==NNTI_GET_SRC) || (remote_op==NNTI_PUT_DST) || (remote_op==(NNTI_GET_SRC|NNTI_PUT_DST))) {
            log_error(nnti_ee_debug_level, "Target buffer ACKs are disabled.  You cannot wait on RDMA target buffers.");
            memset(status, 0, sizeof(NNTI_status_t));
            status->op     = remote_op;
            status->result = NNTI_EINVAL;
            return(NNTI_EINVAL);
        }
    }

    if (is_buf_op_complete(reg_buf) == TRUE) {
        log_debug(nnti_event_debug_level, "buffer op already complete (reg_buf=%p)", reg_buf);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(nnti_event_debug_level, "buffer op NOT complete (reg_buf=%p)", reg_buf);

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(nnti_debug_level, "caught abort signal");
                nnti_rc=NNTI_ECANCELED;
                break;
            }

            check_listen_socket_for_new_connections();

nthread_lock(&nnti_gni_lock);

            if (is_buf_op_complete(reg_buf) == TRUE) {
                log_debug(nnti_event_debug_level, "buffer op already complete (reg_buf=%p)", reg_buf);
                nnti_rc = NNTI_OK;
nthread_unlock(&nnti_gni_lock);
                break;
            }

            cq_hdl=get_cq(reg_buf);

            if ((gni_mem_hdl->type == REQUEST_BUFFER) &&
                (q_hdl->wc_buffer[q_hdl->req_processed].ack_received==1)) {
//                log_debug(nnti_event_debug_level, "processing out of order on cq_hdl(%llu) (l.qw1=%llu l.qw2=%llu r.qw1=%llu r.qw2=%llu)",
//                        (uint64_t)cq_hdl,
//                        (uint64_t)wr->post_desc.local_mem_hndl.qword1,
//                        (uint64_t)wr->post_desc.local_mem_hndl.qword2,
//                        (uint64_t)wr->post_desc.remote_mem_hndl.qword1,
//                        (uint64_t)wr->post_desc.remote_mem_hndl.qword2);
                memset(&ev_data, 0, sizeof(ev_data));
//                memset(&wc, 0, sizeof(nnti_gni_work_completion));
                process_event(reg_buf, cq_hdl, &ev_data);
//                log_debug(nnti_event_debug_level, "out of order processing complete on cq_hdl(%llu) (l.qw1=%llu l.qw2=%llu r.qw1=%llu r.qw2=%llu)",
//                        (uint64_t)cq_hdl,
//                        (uint64_t)wr->post_desc.local_mem_hndl.qword1,
//                        (uint64_t)wr->post_desc.local_mem_hndl.qword2,
//                        (uint64_t)wr->post_desc.remote_mem_hndl.qword1,
//                        (uint64_t)wr->post_desc.remote_mem_hndl.qword2);
                if (is_buf_op_complete(reg_buf) == TRUE) {
                    nnti_rc = NNTI_OK;
nthread_unlock(&nnti_gni_lock);
                    break;
                } else {
nthread_unlock(&nnti_gni_lock);
                    continue;
                }
            } else {
                memset(&ev_data, 0, sizeof(ev_data));
//                log_debug(nnti_event_debug_level, "calling CqWaitEvent(wait) on cq_hdl(%llu) (l.qw1=%llu l.qw2=%llu r.qw1=%llu r.qw2=%llu)",
//                        (uint64_t)cq_hdl,
//                        (uint64_t)wr->post_desc.local_mem_hndl.qword1,
//                        (uint64_t)wr->post_desc.local_mem_hndl.qword2,
//                        (uint64_t)wr->post_desc.remote_mem_hndl.qword1,
//                        (uint64_t)wr->post_desc.remote_mem_hndl.qword2);
//                log_debug(nnti_event_debug_level, "calling CqWaitEvent(wait)");
                trios_start_timer(call_time);
//                nthread_lock(&nnti_gni_lock);
                rc=GNI_CqWaitEvent_wrapper(cq_hdl, timeout_per_call, &ev_data);
//                nthread_unlock(&nnti_gni_lock);
                print_cq_event(&ev_data);
                trios_stop_timer("NNTI_gni_wait - CqWaitEvent", call_time);
                log_debug(nnti_event_debug_level, "CqWaitEvent(wait) complete");
                if (rc!=GNI_RC_SUCCESS) log_debug(nnti_debug_level, "CqWaitEvent() on cq_hdl(%llu) failed: %d", (uint64_t)cq_hdl, rc);
            }

            /* case 1: success */
            if (rc == GNI_RC_SUCCESS) {
                nnti_rc = NNTI_OK;
            }
            /* case 2: timed out */
            else if ((rc==GNI_RC_TIMEOUT) || (rc==GNI_RC_NOT_DONE)) {
                elapsed_time = (trios_get_time_ms() - entry_time);

                /* if the caller asked for a legitimate timeout, we need to exit */
                if (((timeout > 0) && (elapsed_time >= timeout)) || trios_exit_now()) {
                    log_debug(nnti_debug_level, "CqWaitEvent timed out...timeout(%d) elapsed_time(%d) exit_now(%d)",
                            timeout, elapsed_time, trios_exit_now());
                    nnti_rc = NNTI_ETIMEDOUT;
nthread_unlock(&nnti_gni_lock);
                    break;
                }
                /* continue if the timeout has not expired */
                log_debug(nnti_event_debug_level, "CqWaitEvent timedout...retrying...timeout(%d) elapsed_time(%d) exit_now(%d)",
                        timeout, elapsed_time, trios_exit_now());
                //            if (elapsed_time >= 500) {
                //                fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                //                        "NNTI_wait() timeout(5+ sec)", reg_buf);
                //            }

nthread_unlock(&nnti_gni_lock);
                continue;
            }
            /* case 3: failure */
            else {
                char errstr[1024];
                uint32_t recoverable=0;
                GNI_CqErrorStr(ev_data, errstr, 1024);
                GNI_CqErrorRecoverable(ev_data, &recoverable);

                log_error(nnti_debug_level, "CqWaitEvent failed (cq_hdl=%llu ; rc=%d) (recoverable=%llu) : %s",
                        (uint64_t)cq_hdl, rc, (uint64_t)recoverable, errstr);
                print_failed_cq(reg_buf);
                fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                        "NNTI_wait() error", reg_buf);
                nnti_rc = NNTI_EIO;
                break;
            }

            if (rc == GNI_RC_SUCCESS) {
                print_cq_event(&ev_data);

                if (!GNI_CQ_STATUS_OK(ev_data)) {
                    char errstr[1024];
                    GNI_CqErrorStr(ev_data, errstr, 1024);
                    log_error(nnti_debug_level, "Failed status %s (%d)",
                            errstr,
                            GNI_CQ_GET_STATUS(ev_data));
                    nnti_rc=NNTI_EIO;
                    break;
                }

                wait_buf=decode_event_buffer(reg_buf, &ev_data);
//                memset(&wc, 0, sizeof(nnti_gni_work_completion));
                process_event(wait_buf, cq_hdl, &ev_data);
                if (is_buf_op_complete(reg_buf) == TRUE) {
                    nnti_rc = NNTI_OK;
nthread_unlock(&nnti_gni_lock);
                    break;
                }
            }
nthread_unlock(&nnti_gni_lock);
        }
    }

//    print_wc(&wc);
//    if (nnti_rc == NNTI_OK) {
//        print_raw_buf((char *)reg_buf->payload+wr->wc.byte_offset, gni_mem_hdl->last_wc.byte_len);
//    }

    create_status(reg_buf, remote_op, nnti_rc, &ev_data, status);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_wait", status);
    }

    if (nnti_rc==NNTI_OK) {
        switch (gni_mem_hdl->type) {
            case REQUEST_BUFFER:
            case RECEIVE_BUFFER:
                nthread_lock(&gni_mem_hdl->wr_queue_lock);
                wr=gni_mem_hdl->wr_queue.front();
                gni_mem_hdl->wr_queue.pop_front();
                repost_recv_work_request((NNTI_buffer_t *)reg_buf, wr);
                nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                break;
            case GET_SRC_BUFFER:
            case PUT_DST_BUFFER:
            case RDMA_TARGET_BUFFER:
                if (config.use_rdma_target_ack) {
                    nthread_lock(&gni_mem_hdl->wr_queue_lock);
                    wr=gni_mem_hdl->wr_queue.front();
                    gni_mem_hdl->wr_queue.pop_front();
                    repost_recv_work_request((NNTI_buffer_t *)reg_buf, wr);
                    nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                }
                break;
            case SEND_BUFFER:
            case GET_DST_BUFFER:
            case PUT_SRC_BUFFER:
                nthread_lock(&gni_mem_hdl->wr_queue_lock);
                wr=gni_mem_hdl->wr_queue.front();
                gni_mem_hdl->wr_queue.pop_front();
                nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                del_wr_wrhash(wr);
                if (config.use_wr_pool) {
                    wr_pool_initiator_push(wr);
                } else {
                    unregister_wc(wr);
                    free(wr);
                }
                break;
            case UNKNOWN_BUFFER:
            default:
                log_error(nnti_debug_level, "unknown buffer type(%llu).", gni_mem_hdl->type);
                break;
        }
    }

    log_debug(nnti_ee_debug_level, "exit");
    return((NNTI_result_t)nnti_rc);
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
NNTI_result_t NNTI_gni_waitany (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status)
{
    int nnti_rc=NNTI_OK;
    gni_memory_handle        *gni_mem_hdl=NULL;
    gni_work_request         *wr=NULL;
    gni_request_queue_handle *q_hdl=NULL;
    gni_connection           *conn=NULL;
    const NNTI_buffer_t  *wait_buf=NULL;

    uint32_t which_cq=0;

//    gni_cq_handle_t cq_hdl=transport_global_data.ep_cq_hdl;
    gni_cq_entry_t  ev_data;

    nnti_gni_work_completion wc;

    gni_return_t rc=GNI_RC_SUCCESS;
    long elapsed_time = 0;
    long timeout_per_call;

    long entry_time=trios_get_time_ms();

    trios_declare_timer(call_time);

    log_debug(nnti_ee_debug_level, "enter");

    check_listen_socket_for_new_connections();

    if (!config.use_rdma_target_ack) {
        if ((remote_op==NNTI_GET_SRC) || (remote_op==NNTI_PUT_DST) || (remote_op==(NNTI_GET_SRC|NNTI_PUT_DST))) {
            log_error(nnti_ee_debug_level, "Target buffer ACKs are disabled.  You cannot wait on RDMA target buffers.");
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
                assert(((gni_memory_handle *)buf_list[i]->transport_private)->type != REQUEST_BUFFER);
            }
        }
    }
    assert(status);

    if (buf_count == 1) {
        nnti_rc=NNTI_gni_wait(buf_list[0], remote_op, timeout, status);
        *which=0;
        goto cleanup;
    }

    if (is_any_buf_op_complete(buf_list, buf_count, which) == TRUE) {
        log_debug(nnti_event_debug_level, "buffer op already complete (which=%u, buf_list[%d]=%p)", *which, *which, buf_list[*which]);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(nnti_event_debug_level, "buffer op NOT complete (buf_list=%p)", buf_list);

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(nnti_debug_level, "caught abort signal");
                nnti_rc=NNTI_ECANCELED;
                break;
            }

            check_listen_socket_for_new_connections();

            if (is_any_buf_op_complete(buf_list, buf_count, which) == TRUE) {
                log_debug(nnti_event_debug_level, "buffer op already complete (which=%u, buf_list[%d]=%p)", *which, *which, buf_list[*which]);
                nnti_rc = NNTI_OK;
                break;
            }

            memset(&ev_data, 0, sizeof(ev_data));
            log_debug(nnti_event_debug_level, "calling CqVectorWaitEvent(wait)");
//            log_debug(nnti_event_debug_level, "calling CqWaitEvent(wait) on cq_hdl(%llu) (l.qw1=%llu l.qw2=%llu r.qw1=%llu r.qw2=%llu)",
//                    (uint64_t)cq_hdl,
//                    (uint64_t)wr->post_desc.local_mem_hndl.qword1,
//                    (uint64_t)wr->post_desc.local_mem_hndl.qword2,
//                    (uint64_t)wr->post_desc.remote_mem_hndl.qword1,
//                    (uint64_t)wr->post_desc.remote_mem_hndl.qword2);
            //            log_debug(nnti_event_debug_level, "calling CqWaitEvent(wait)");
            trios_start_timer(call_time);
            //            nthread_lock(&nnti_gni_lock);
            rc=GNI_CqVectorWaitEvent (&transport_global_data.cq_list[0], 2, timeout_per_call, &ev_data, &which_cq);
            print_cq_event(&ev_data);
            //            nthread_unlock(&nnti_gni_lock);
            trios_stop_timer("NNTI_gni_waitany - CqVectorWaitEvent", call_time);
            log_debug(nnti_event_debug_level, "CqVectorWaitEvent(wait) complete - which_cq=%lu, cq_list[%lu]=%llu",
                    (unsigned long)which_cq, (unsigned long)which_cq, (unsigned long long)transport_global_data.cq_list[which_cq]);
            if (rc!=GNI_RC_SUCCESS) log_debug(nnti_debug_level, "CqVectorWaitEvent() failed: %d", rc);

            /* case 1: success */
            if (rc == GNI_RC_SUCCESS) {
                nnti_rc = NNTI_OK;
            }
            /* case 2: timed out */
            else if ((rc==GNI_RC_TIMEOUT) || (rc==GNI_RC_NOT_DONE)) {
                elapsed_time = (trios_get_time_ms() - entry_time);

                /* if the caller asked for a legitimate timeout, we need to exit */
                if (((timeout > 0) && (elapsed_time >= timeout)) || trios_exit_now()) {
                    log_debug(nnti_debug_level, "CqVectorWaitEvent timed out...timeout(%d) elapsed_time(%d) exit_now(%d)",
                            timeout, elapsed_time, trios_exit_now());
                    nnti_rc = NNTI_ETIMEDOUT;
                    break;
                }
                /* continue if the timeout has not expired */
                log_debug(nnti_event_debug_level, "CqVectorWaitEvent timedout...retrying...timeout(%d) elapsed_time(%d) exit_now(%d)",
                        timeout, elapsed_time, trios_exit_now());
                //            if (elapsed_time >= 500) {
                //                fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                //                        "NNTI_wait() timeout(5+ sec)", reg_buf);
                //            }

                continue;
            }
            /* case 3: failure */
            else {
                char errstr[1024];
                uint32_t recoverable=0;
                GNI_CqErrorStr(ev_data, errstr, 1024);
                GNI_CqErrorRecoverable(ev_data, &recoverable);

                log_error(nnti_debug_level, "CqVectorWaitEvent failed (rc=%d) (recoverable=%llu) : %s",
                        rc, (uint64_t)recoverable, errstr);
//                print_failed_cq(reg_buf);
//                fprint_NNTI_buffer(logger_get_file(), "reg_buf",
//                        "NNTI_wait() error", reg_buf);
                nnti_rc = NNTI_EIO;
                break;
            }

            if (rc == GNI_RC_SUCCESS) {
                print_cq_event(&ev_data);

                if (!GNI_CQ_STATUS_OK(ev_data)) {
                    char errstr[1024];
                    GNI_CqErrorStr(ev_data, errstr, 1024);
                    log_error(nnti_debug_level, "Failed status %s (%d)",
                            errstr,
                            GNI_CQ_GET_STATUS(ev_data));
                    nnti_rc=NNTI_EIO;
                    break;
                }

                memset(&wc, 0, sizeof(nnti_gni_work_completion));
                wait_buf=decode_event_buffer(buf_list[0], &ev_data);
//                process_event(wait_buf, transport_global_data.cq_list[which_cq], &ev_data, &wc);
                process_event(wait_buf, transport_global_data.cq_list[which_cq], &ev_data);
                if (is_any_buf_op_complete(buf_list, buf_count, which) == TRUE) {
                    nnti_rc = NNTI_OK;
                    break;
                }
            }
        }
    }

    create_status(buf_list[*which], remote_op, nnti_rc, NULL, status);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_waitany", status);
    }

    if (nnti_rc==NNTI_OK) {
        gni_mem_hdl=(gni_memory_handle *)buf_list[*which]->transport_private;
        assert(gni_mem_hdl);
        switch (gni_mem_hdl->type) {
            case REQUEST_BUFFER:
            case RECEIVE_BUFFER:
                nthread_lock(&gni_mem_hdl->wr_queue_lock);
                wr=gni_mem_hdl->wr_queue.front();
                gni_mem_hdl->wr_queue.pop_front();
                repost_recv_work_request((NNTI_buffer_t *)buf_list[*which], wr);
                nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                break;
            case GET_SRC_BUFFER:
            case PUT_DST_BUFFER:
            case RDMA_TARGET_BUFFER:
                if (config.use_rdma_target_ack) {
                    nthread_lock(&gni_mem_hdl->wr_queue_lock);
                    wr=gni_mem_hdl->wr_queue.front();
                    gni_mem_hdl->wr_queue.pop_front();
                    repost_recv_work_request((NNTI_buffer_t *)buf_list[*which], wr);
                    nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                }
                break;
            case SEND_BUFFER:
            case GET_DST_BUFFER:
            case PUT_SRC_BUFFER:
                nthread_lock(&gni_mem_hdl->wr_queue_lock);
                wr=gni_mem_hdl->wr_queue.front();
                gni_mem_hdl->wr_queue.pop_front();
                nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                del_wr_wrhash(wr);
                if (config.use_wr_pool) {
                    wr_pool_initiator_push(wr);
                } else {
                    unregister_wc(wr);
                    free(wr);
                }
                break;
            case UNKNOWN_BUFFER:
            default:
                log_error(nnti_debug_level, "unknown buffer type(%llu).", gni_mem_hdl->type);
                break;
        }
    }

cleanup:
    log_debug(nnti_ee_debug_level, "exit");

    return((NNTI_result_t)nnti_rc);
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
NNTI_result_t NNTI_gni_waitall (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t       **status)
{
    int nnti_rc=NNTI_OK;
    gni_memory_handle        *gni_mem_hdl=NULL;
    gni_work_request         *wr=NULL;
    gni_request_queue_handle *q_hdl=NULL;
    gni_connection           *conn=NULL;
    const NNTI_buffer_t  *wait_buf=NULL;

    uint32_t which_cq=0;

//    gni_cq_handle_t cq_hdl=transport_global_data.ep_cq_hdl;
    gni_cq_entry_t  ev_data;

    nnti_gni_work_completion wc;

    gni_return_t rc=GNI_RC_SUCCESS;
    long elapsed_time = 0;
    long timeout_per_call;

    long entry_time=trios_get_time_ms();

    trios_declare_timer(call_time);

    log_debug(nnti_ee_debug_level, "enter");

    check_listen_socket_for_new_connections();

    if (!config.use_rdma_target_ack) {
        if ((remote_op==NNTI_GET_SRC) || (remote_op==NNTI_PUT_DST) || (remote_op==(NNTI_GET_SRC|NNTI_PUT_DST))) {
            log_error(nnti_ee_debug_level, "Target buffer ACKs are disabled.  You cannot wait on RDMA target buffers.");
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
                assert(((gni_memory_handle *)buf_list[i]->transport_private)->type != REQUEST_BUFFER);
            }
        }
    }
    assert(status);

    if (buf_count == 1) {
        nnti_rc=NNTI_gni_wait(buf_list[0], remote_op, timeout, status[0]);
        goto cleanup;
    }

    if (is_all_buf_ops_complete(buf_list, buf_count) == TRUE) {
        log_debug(nnti_event_debug_level, "all buffer ops already complete");
        nnti_rc = NNTI_OK;
    } else {
        log_debug(nnti_event_debug_level, "all buffer ops NOT complete");
//        reset_op_state(reg_buf);

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(nnti_debug_level, "caught abort signal");
                nnti_rc=NNTI_ECANCELED;
                break;
            }

            check_listen_socket_for_new_connections();

            if (is_all_buf_ops_complete(buf_list, buf_count) == TRUE) {
                log_debug(nnti_event_debug_level, "all buffer ops already complete");
                nnti_rc = NNTI_OK;
                break;
            }

            memset(&ev_data, 0, sizeof(ev_data));
            log_debug(nnti_event_debug_level, "calling CqVectorWaitEvent(wait)");
//            log_debug(nnti_event_debug_level, "calling CqWaitEvent(wait) on cq_hdl(%llu) (l.qw1=%llu l.qw2=%llu r.qw1=%llu r.qw2=%llu)",
//                    (uint64_t)cq_hdl,
//                    (uint64_t)wr->post_desc.local_mem_hndl.qword1,
//                    (uint64_t)wr->post_desc.local_mem_hndl.qword2,
//                    (uint64_t)wr->post_desc.remote_mem_hndl.qword1,
//                    (uint64_t)wr->post_desc.remote_mem_hndl.qword2);
            //            log_debug(nnti_event_debug_level, "calling CqWaitEvent(wait)");
            trios_start_timer(call_time);
            //            nthread_lock(&nnti_gni_lock);
            rc=GNI_CqVectorWaitEvent (&transport_global_data.cq_list[0], 2, timeout_per_call, &ev_data, &which_cq);
            print_cq_event(&ev_data);
            //            nthread_unlock(&nnti_gni_lock);
            trios_stop_timer("NNTI_gni_waitall - CqVectorWaitEvent", call_time);
            log_debug(nnti_event_debug_level, "CqVectorWaitEvent(wait) complete - which_cq=%lu, cq_list[%lu]=%llu",
                    (unsigned long)which_cq, (unsigned long)which_cq, (unsigned long long)transport_global_data.cq_list[which_cq]);
            if (rc!=GNI_RC_SUCCESS) log_debug(nnti_debug_level, "CqVectorWaitEvent() failed: %d", rc);

            /* case 1: success */
            if (rc == GNI_RC_SUCCESS) {
                nnti_rc = NNTI_OK;
            }
            /* case 2: timed out */
            else if ((rc==GNI_RC_TIMEOUT) || (rc==GNI_RC_NOT_DONE)) {
                elapsed_time = (trios_get_time_ms() - entry_time);

                /* if the caller asked for a legitimate timeout, we need to exit */
                if (((timeout > 0) && (elapsed_time >= timeout)) || trios_exit_now()) {
                    log_debug(nnti_debug_level, "CqVectorWaitEvent timed out...timeout(%d) elapsed_time(%d) exit_now(%d)",
                            timeout, elapsed_time, trios_exit_now());
                    nnti_rc = NNTI_ETIMEDOUT;
                    break;
                }
                /* continue if the timeout has not expired */
                log_debug(nnti_event_debug_level, "CqVectorWaitEvent timedout...retrying...timeout(%d) elapsed_time(%d) exit_now(%d)",
                        timeout, elapsed_time, trios_exit_now());
                //            if (elapsed_time >= 500) {
                //                fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                //                        "NNTI_wait() timeout(5+ sec)", reg_buf);
                //            }

                continue;
            }
            /* case 3: failure */
            else {
                char errstr[1024];
                uint32_t recoverable=0;
                GNI_CqErrorStr(ev_data, errstr, 1024);
                GNI_CqErrorRecoverable(ev_data, &recoverable);

                log_error(nnti_debug_level, "CqVectorWaitEvent failed (rc=%d) (recoverable=%llu) : %s",
                        rc, (uint64_t)recoverable, errstr);
//                print_failed_cq(reg_buf);
//                fprint_NNTI_buffer(logger_get_file(), "reg_buf",
//                        "NNTI_wait() error", reg_buf);
                nnti_rc = NNTI_EIO;
                break;
            }

            if (rc == GNI_RC_SUCCESS) {
                print_cq_event(&ev_data);

                if (!GNI_CQ_STATUS_OK(ev_data)) {
                    char errstr[1024];
                    GNI_CqErrorStr(ev_data, errstr, 1024);
                    log_error(nnti_debug_level, "Failed status %s (%d)",
                            errstr,
                            GNI_CQ_GET_STATUS(ev_data));
                    nnti_rc=NNTI_EIO;
                    break;
                }

                memset(&wc, 0, sizeof(nnti_gni_work_completion));
                wait_buf=decode_event_buffer(buf_list[0], &ev_data);
//                process_event(wait_buf, transport_global_data.cq_list[which_cq], &ev_data, &wc);
                process_event(wait_buf, transport_global_data.cq_list[which_cq], &ev_data);
                if (is_all_buf_ops_complete(buf_list, buf_count) == TRUE) {
                    nnti_rc = NNTI_OK;
                    break;
                }
            }
        }
    }

    for (int i=0;i<buf_count;i++) {
        create_status(buf_list[i], remote_op, nnti_rc, NULL, status[i]);

        if (logging_debug(nnti_debug_level)) {
            fprint_NNTI_status(logger_get_file(), "status",
                    "end of NNTI_waitall", status[i]);
        }

        if (nnti_rc==NNTI_OK) {
            gni_mem_hdl=(gni_memory_handle *)buf_list[i]->transport_private;
            assert(gni_mem_hdl);
            switch (gni_mem_hdl->type) {
                case REQUEST_BUFFER:
                case RECEIVE_BUFFER:
                    nthread_lock(&gni_mem_hdl->wr_queue_lock);
                    wr=gni_mem_hdl->wr_queue.front();
                    gni_mem_hdl->wr_queue.pop_front();
                    repost_recv_work_request((NNTI_buffer_t *)buf_list[i], wr);
                    nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                    break;
                case GET_SRC_BUFFER:
                case PUT_DST_BUFFER:
                case RDMA_TARGET_BUFFER:
                    if (config.use_rdma_target_ack) {
                        nthread_lock(&gni_mem_hdl->wr_queue_lock);
                        wr=gni_mem_hdl->wr_queue.front();
                        gni_mem_hdl->wr_queue.pop_front();
                        repost_recv_work_request((NNTI_buffer_t *)buf_list[i], wr);
                        nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                    }
                    break;
                case SEND_BUFFER:
                case GET_DST_BUFFER:
                case PUT_SRC_BUFFER:
                    nthread_lock(&gni_mem_hdl->wr_queue_lock);
                    wr=gni_mem_hdl->wr_queue.front();
                    gni_mem_hdl->wr_queue.pop_front();
                    nthread_unlock(&gni_mem_hdl->wr_queue_lock);
                    del_wr_wrhash(wr);
                    if (config.use_wr_pool) {
                        wr_pool_initiator_push(wr);
                    } else {
                        unregister_wc(wr);
                        free(wr);
                    }
                    break;
                case UNKNOWN_BUFFER:
                default:
                    log_error(nnti_debug_level, "unknown buffer type(%llu).", gni_mem_hdl->type);
                    break;
            }
        }
    }

cleanup:
    log_debug(nnti_ee_debug_level, "exit");

    return((NNTI_result_t)nnti_rc);
}

/**
 * @brief Disable this transport.
 *
 * Shutdown the transport.  Any outstanding sends, gets and puts will be
 * canceled.  Any new transport requests will fail.
 *
 */
NNTI_result_t NNTI_gni_fini (
        const NNTI_transport_t *trans_hdl)
{
    int rc=GNI_RC_SUCCESS; /* return code */

    log_debug(nnti_ee_debug_level, "enter");

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

    rc=GNI_CdmDestroy(transport_global_data.cdm_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "CdmCreate() failed: %d", rc);
        rc=NNTI_EINVAL;
    }

    nthread_lock_fini(&nnti_gni_lock);
    nthread_lock_fini(&nnti_conn_peer_lock);
    nthread_lock_fini(&nnti_conn_instance_lock);
    nthread_lock_fini(&nnti_wr_wrhash_lock);
    nthread_lock_fini(&nnti_buf_bufhash_lock);
    nthread_lock_fini(&nnti_wr_pool_lock);

    gni_initialized = false;

    log_debug(nnti_ee_debug_level, "exit");

    return(NNTI_OK);
}

static NNTI_result_t register_wc(gni_work_request *wr)
{
    int rc=GNI_RC_SUCCESS; /* return code */

    int cq_cnt=1;

    trios_declare_timer(call_time);

    gni_connection *conn=NULL;
    gni_memory_handle *hdl=NULL;
    gni_cq_handle_t wc_mem_cq_hdl;

    log_debug(nnti_debug_level, "enter wr(%p)", wr);

    assert(wr);

    hdl=(gni_memory_handle *)wr->reg_buf->transport_private;
    assert(hdl);

    wc_mem_cq_hdl=NULL;

    if (need_wc_mem_cq(hdl) == 1) {
        wc_mem_cq_hdl=transport_global_data.mem_cq_hdl;
    }

    trios_start_timer(call_time);
    rc=GNI_MemRegister (transport_global_data.nic_hdl,
            (uint64_t)&wr->wc,
            sizeof(nnti_gni_work_completion),
            wc_mem_cq_hdl,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &wr->wc_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "MemRegister(wc_mem_hdl) failed: rc=%d, %s", rc, strerror(errno));
        goto cleanup;
    }
    trios_stop_timer("wc register", call_time);

    wr->wc_registered=TRUE;

    log_debug(nnti_debug_level, "register wc_mem_cq_hdl  =%llu", (uint64_t)wc_mem_cq_hdl);
    log_debug(nnti_debug_level, "register hdl->wc_mem_hdl=(%llu,%llu)", (uint64_t)wr->wc_mem_hdl.qword1, (uint64_t)wr->wc_mem_hdl.qword2);

    log_debug(nnti_debug_level, "exit  wr(%p)", wr);

    return((NNTI_result_t)GNI_RC_SUCCESS);

cleanup:
//    unregister_memory(hdl);

    switch(rc) {
        case GNI_RC_SUCCESS:
            rc=(int)NNTI_OK;
        default:
            rc=(int)NNTI_EIO;
    }

    return ((NNTI_result_t)rc);
}

static NNTI_result_t register_memory(gni_memory_handle *hdl, void *buf, uint64_t len)
{
    int rc=GNI_RC_SUCCESS; /* return code */

    int cq_cnt=1;

    uint32_t flags=GNI_MEM_READWRITE;

    trios_declare_timer(call_time);

    gni_connection *conn=NULL;

    gni_cq_handle_t mem_cq_hdl;
    gni_cq_handle_t wc_mem_cq_hdl;

    assert(hdl);

    mem_cq_hdl   =NULL;
    wc_mem_cq_hdl=NULL;

    log_debug(nnti_debug_level, "enter hdl(%p) buffer(%p) len(%d)", hdl, buf, len);

    if (config.pi_ordering==PI_ORDERING_STRICT) {
        log_debug(nnti_debug_level, "using STRICT ordering");
        flags = GNI_MEM_READWRITE|GNI_MEM_STRICT_PI_ORDERING;
    } else if (config.pi_ordering==PI_ORDERING_RELAXED) {
        log_debug(nnti_debug_level, "using RELAXED ordering");
        flags = GNI_MEM_READWRITE|GNI_MEM_RELAXED_PI_ORDERING;
    } else if (config.pi_ordering==PI_ORDERING_DEFAULT) {
        log_debug(nnti_debug_level, "using DEFAULT ordering");
        flags = GNI_MEM_READWRITE;
    }

    if (need_mem_cq(hdl) == 1) {
        mem_cq_hdl=transport_global_data.mem_cq_hdl;
    }

    trios_start_timer(call_time);
    rc=GNI_MemRegister (transport_global_data.nic_hdl,
            (uint64_t)buf,
            len,
            mem_cq_hdl,
            flags,
            (uint32_t)-1,
            &hdl->mem_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "MemRegister(mem_hdl) failed: rc=%d, %s", rc, strerror(errno));
        goto cleanup;
    }
    trios_stop_timer("buf register", call_time);

    log_debug(nnti_debug_level, "register mem_cq_hdl     =%llu", (uint64_t)mem_cq_hdl);
    log_debug(nnti_debug_level, "register hdl->mem_hdl   =(%llu,%llu)", (uint64_t)hdl->mem_hdl.qword1, (uint64_t)hdl->mem_hdl.qword2);

    log_debug(nnti_debug_level, "exit  hdl(%p) buffer(%p)", hdl, buf);

    return((NNTI_result_t)GNI_RC_SUCCESS);

cleanup:
//    unregister_memory(hdl);

    switch(rc) {
        case GNI_RC_SUCCESS:
            rc=(int)NNTI_OK;
        default:
            rc=(int)NNTI_EIO;
    }

    return ((NNTI_result_t)rc);
}

static NNTI_result_t unregister_wc(gni_work_request *wr)
{
    int rc=GNI_RC_SUCCESS; /* return code */

    int cq_cnt=1;

    trios_declare_timer(call_time);

    gni_connection *conn=NULL;
    gni_memory_handle *hdl=NULL;
    gni_cq_handle_t wc_mem_cq_hdl;

    log_debug(nnti_debug_level, "enter wr(%p)", wr);

    assert(wr);

    if (wr->wc_registered==FALSE) {
        log_debug(nnti_debug_level, "exit wr(%p) - not registered", wr);
        return(NNTI_OK);
    }

    hdl=(gni_memory_handle *)wr->reg_buf->transport_private;
    assert(hdl);

    log_debug(nnti_debug_level, "unregister hdl->wc_mem_hdl   =(%llu,%llu)", (uint64_t)wr->wc_mem_hdl.qword1, (uint64_t)wr->wc_mem_hdl.qword2);

    trios_start_timer(call_time);
    rc=GNI_MemDeregister (transport_global_data.nic_hdl, &wr->wc_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "MemDeregister(wc_mem_hdl) failed: %d", rc);
    }
    trios_stop_timer("wc deregister", call_time);

    switch(rc) {
        case GNI_RC_SUCCESS:
            rc=NNTI_OK;
        default:
            rc=NNTI_EIO;
    }

    wr->wc_registered=FALSE;

    log_debug(nnti_debug_level, "exit wr(%p)", wr);

    return ((NNTI_result_t)rc);
}

static NNTI_result_t unregister_memory(gni_memory_handle *hdl)
{
    int rc=GNI_RC_SUCCESS; /* return code */
    trios_declare_timer(call_time);

    gni_connection *conn=NULL;

    assert(hdl);

    log_debug(nnti_debug_level, "enter hdl(%p)", hdl);

    log_debug(nnti_debug_level, "unregister hdl->mem_hdl   =(%llu,%llu)", (uint64_t)hdl->mem_hdl.qword1, (uint64_t)hdl->mem_hdl.qword2);

    trios_start_timer(call_time);
    rc=GNI_MemDeregister (transport_global_data.nic_hdl, &hdl->mem_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "MemDeregister(mem_hdl) failed: %d", rc);
    }
    trios_stop_timer("buf deregister", call_time);

    switch(rc) {
        case GNI_RC_SUCCESS:
            rc=NNTI_OK;
        default:
            rc=NNTI_EIO;
    }

    log_debug(nnti_debug_level, "exit  hdl(%p)", hdl);

    return ((NNTI_result_t)rc);
}

//if (config.use_rdma_target_ack) {
static void send_rdma_wc (
        gni_work_request    *wr,
        const NNTI_buffer_t *local_buf,
        const NNTI_buffer_t *remote_buf)
{
    int rc=NNTI_OK;
    trios_declare_timer(call_time);

    gni_ep_handle_t        ep_hdl;

    gni_memory_handle *gni_mem_hdl=NULL;

    bool use_fma=true;

    assert(local_buf);
    assert(remote_buf);

    log_debug(nnti_ee_debug_level, "enter");

    gni_mem_hdl=(gni_memory_handle *)local_buf->transport_private;
    assert(gni_mem_hdl);

    print_wc(&wr->wc);

    memset(&wr->wc_post_desc, 0, sizeof(gni_post_descriptor_t));
    set_wc_post_desc(&wr->wc_post_desc, sizeof(nnti_gni_work_completion));

    wr->wc_post_desc.local_addr     =(uint64_t)&wr->wc;
    wr->wc_post_desc.local_mem_hndl =wr->wc_mem_hdl;
    wr->wc_post_desc.remote_addr    =wr->wc_dest_addr;
    wr->wc_post_desc.remote_mem_hndl=wr->wc_dest_mem_hdl;
    wr->wc_post_desc.length         =sizeof(nnti_gni_work_completion);

    gni_connection *conn=get_conn_peer(&remote_buf->buffer_owner);
    assert(conn);
    ep_hdl=conn->ep_hdl;

    log_debug(nnti_debug_level, "sending ACK to (%s, ep_hdl=%llu, peer_instance=%lu)",
            remote_buf->buffer_owner.url, ep_hdl, remote_buf->buffer_owner.peer.NNTI_remote_process_t_u.gni.inst_id);

    if (config.rdma_mode==RDMA_CROSSOVER) {
        if (wr->wc_post_desc.length > config.fma_bte_crossover_size) {
            use_fma=false;
        } else {
            use_fma=true;
        }
    } else if (config.rdma_mode==RDMA_BTE) {
        use_fma=false;
    } else if ((config.rdma_mode==RDMA_FMA) || (config.rdma_mode==RDMA_MIXED)) {
        use_fma=true;
    }

    if (use_fma) {
        log_debug(nnti_debug_level, "calling PostFma(fma wc ep_hdl(%llu) transport_global_data.cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                ep_hdl, transport_global_data.ep_cq_hdl, wr->wc_post_desc.local_addr, wr->wc_post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostFma(ep_hdl, &wr->wc_post_desc);
        trios_stop_timer("PostFma ack", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostFma(fma wc) failed: %d", rc);
    } else {
        log_debug(nnti_debug_level, "calling PostRdma(rdma wc ep_hdl(%llu) transport_global_data.cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                ep_hdl, transport_global_data.ep_cq_hdl, wr->wc_post_desc.local_addr, wr->wc_post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostRdma(ep_hdl, &wr->wc_post_desc);
        trios_stop_timer("PostRdma ack", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostRdma(rdma wc) failed: %d", rc);
    }

    log_debug(nnti_ee_debug_level, "exit");

    return;
}
//#endif


static int need_mem_cq(gni_memory_handle *gni_mem_hdl)
{
    int need_cq=0;

    gni_work_request *wr=NULL;

    log_debug(nnti_ee_debug_level, "enter");

    assert(gni_mem_hdl);

    if (!config.use_rdma_events) {
        need_cq=0;
    } else {
        switch (gni_mem_hdl->type) {
        case GET_DST_BUFFER:
        case PUT_SRC_BUFFER:
        case SEND_BUFFER:
            need_cq=0;
            break;
        case GET_SRC_BUFFER:
            if (config.rdma_mode==RDMA_CROSSOVER) {
                need_cq=0;
            } else {
                need_cq=1;
            }
            break;
        case PUT_DST_BUFFER:
        case RECEIVE_BUFFER:
        case REQUEST_BUFFER:
            need_cq=1;
            break;
        case RDMA_TARGET_BUFFER:
            if (config.rdma_mode==RDMA_CROSSOVER) {
                need_cq=0;
            } else {
                need_cq=1;
            }
            break;
        case UNKNOWN_BUFFER:
        default:
            need_cq=0;
            break;
        }
    }

    log_debug(nnti_ee_debug_level, "exit(%d)", need_cq);

    return(need_cq);
}

static int need_wc_mem_cq(gni_memory_handle *gni_mem_hdl)
{
    int need_cq=0;

    gni_work_request *wr=NULL;

    log_debug(nnti_ee_debug_level, "enter");

    assert(gni_mem_hdl);

    if (!config.use_rdma_target_ack) {
        switch (gni_mem_hdl->type) {
            case RDMA_TARGET_BUFFER:
            case GET_SRC_BUFFER:
            case PUT_DST_BUFFER:
            case GET_DST_BUFFER:
            case PUT_SRC_BUFFER:
            case SEND_BUFFER:
                need_cq=0;
                break;
            case RECEIVE_BUFFER:
            case REQUEST_BUFFER:
                need_cq=1;
                break;
            case UNKNOWN_BUFFER:
            default:
                need_cq=0;
                break;
        }
    } else {
        switch (gni_mem_hdl->type) {
            case GET_DST_BUFFER:
            case PUT_SRC_BUFFER:
            case SEND_BUFFER:
                need_cq=0;
                break;
            case RDMA_TARGET_BUFFER:
            case GET_SRC_BUFFER:
            case PUT_DST_BUFFER:
            case RECEIVE_BUFFER:
            case REQUEST_BUFFER:
                need_cq=1;
                break;
            case UNKNOWN_BUFFER:
            default:
                need_cq=0;
                break;
        }
    }

    log_debug(nnti_ee_debug_level, "exit(%d)", need_cq);

    return(need_cq);
}

static gni_cq_handle_t get_cq(const NNTI_buffer_t *reg_buf)
{
    gni_memory_handle *gni_mem_hdl=NULL;
    gni_work_request  *wr=NULL;
    gni_connection    *conn=NULL;

    gni_cq_handle_t           cq_hdl=0;
    gni_request_queue_handle *q_hdl=NULL;

    log_debug(nnti_ee_debug_level, "enter");

    assert(reg_buf);

    q_hdl      =&transport_global_data.req_queue;
    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    wr=first_incomplete_wr(gni_mem_hdl);
    assert(wr);

    switch (gni_mem_hdl->type) {
        case GET_DST_BUFFER:
        case PUT_SRC_BUFFER:
            cq_hdl=transport_global_data.ep_cq_hdl;
            break;
        case SEND_BUFFER:
            if (wr->last_op==GNI_OP_SEND_REQUEST) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    conn=get_conn_instance(wr->peer_instance);
                    assert(conn);
                    cq_hdl=conn->queue_local_attrs.req_cq_hdl;
                } else {
                    cq_hdl=transport_global_data.ep_cq_hdl;
                }
            } else if (wr->last_op==GNI_OP_SEND_BUFFER) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    cq_hdl=transport_global_data.ep_cq_hdl;
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    cq_hdl=transport_global_data.ep_cq_hdl;
                }
            } else if (wr->last_op==GNI_OP_PUT_INITIATOR) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    cq_hdl=transport_global_data.ep_cq_hdl;
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    cq_hdl=transport_global_data.ep_cq_hdl;
                }
            }
            break;
        case GET_SRC_BUFFER:
        case PUT_DST_BUFFER:
        case RECEIVE_BUFFER:
            cq_hdl=transport_global_data.mem_cq_hdl;
            break;
        case RDMA_TARGET_BUFFER:
            if ((wr->last_op==GNI_OP_GET_INITIATOR) ||
                (wr->last_op==GNI_OP_PUT_INITIATOR)) {

                cq_hdl=transport_global_data.ep_cq_hdl;
            } else {
                cq_hdl=transport_global_data.mem_cq_hdl;
            }
            break;
        case REQUEST_BUFFER:
            cq_hdl=q_hdl->wc_mem_cq_hdl;
            break;
        case UNKNOWN_BUFFER:
        default:
            log_error(nnti_debug_level, "unknown buffer type(%llu).", gni_mem_hdl->type);
            cq_hdl=(gni_cq_handle_t)-1;
            break;
    }

    switch (gni_mem_hdl->type) {
        case PUT_SRC_BUFFER:
            if (wr->op_state==RDMA_WRITE_INIT) {
                log_debug(nnti_cq_debug_level, "rdma initiator.  waiting on ep_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
            } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                log_debug(nnti_cq_debug_level, "rdma initiator.  waiting on wc_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
            }
            break;
        case GET_DST_BUFFER:
            if (wr->op_state==RDMA_READ_INIT) {
                log_debug(nnti_cq_debug_level, "rdma initiator.  waiting on ep_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
            } else if (wr->op_state==RDMA_READ_NEED_ACK) {
                log_debug(nnti_cq_debug_level, "rdma initiator.  waiting on wc_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
            }
            break;
        case SEND_BUFFER:
            if (wr->last_op==GNI_OP_SEND_REQUEST) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    conn=get_conn_instance(wr->peer_instance);
                    assert(conn);
                    log_debug(nnti_cq_debug_level, "send request using send op.  waiting on req_cq_hdl(%llu).", (uint64_t)conn->queue_local_attrs.req_cq_hdl);
                } else {
                    log_debug(nnti_cq_debug_level, "send request using send op.  waiting on wc_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
                }
            } else if (wr->last_op==GNI_OP_SEND_BUFFER) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(nnti_cq_debug_level, "send buffer using send op.  waiting on ep_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    log_debug(nnti_cq_debug_level, "send buffer using send op.  waiting on wc_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
                }
            } else if (wr->last_op==GNI_OP_PUT_INITIATOR) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(nnti_cq_debug_level, "send buffer using put op.  waiting on ep_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    log_debug(nnti_cq_debug_level, "send buffer using put op.  waiting on wc_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
                }
            }
            break;
        case GET_SRC_BUFFER:
            if (config.use_rdma_events) {
                if (wr->op_state==RDMA_READ_INIT) {
                    log_debug(nnti_cq_debug_level, "rdma target (get_src).  waiting on mem_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
                } else if (wr->op_state==RDMA_READ_NEED_ACK) {
                    log_debug(nnti_cq_debug_level, "rdma target (get_src).  waiting on wc_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
                }
            } else {
                log_debug(nnti_cq_debug_level, "rdma target (get_src).  waiting on wc_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
            }
            break;
        case PUT_DST_BUFFER:
            if (config.use_rdma_events) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(nnti_cq_debug_level, "rdma target (put_dest).  waiting on mem_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    log_debug(nnti_cq_debug_level, "rdma target (put_dest).  waiting on wc_mem_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
                }
            } else {
                log_debug(nnti_cq_debug_level, "rdma target (put_dest/result/receive).  waiting on wc_mem_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
            }
            break;
        case RECEIVE_BUFFER:
            if (config.use_rdma_events) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(nnti_cq_debug_level, "rdma target (result/receive).  waiting on mem_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    log_debug(nnti_cq_debug_level, "rdma target (result/receive).  waiting on wc_mem_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
                }
            } else {
                log_debug(nnti_cq_debug_level, "rdma target (result/receive).  waiting on wc_mem_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
            }
            break;
        case RDMA_TARGET_BUFFER:
            if ((wr->last_op==GNI_OP_GET_INITIATOR) ||
                (wr->last_op==GNI_OP_PUT_INITIATOR)) {

                if (wr->op_state==RDMA_TARGET_INIT) {
                    log_debug(nnti_cq_debug_level, "rdma target (generic target).  waiting on ep_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
                } else if (wr->op_state==RDMA_TARGET_NEED_ACK) {
                    log_debug(nnti_cq_debug_level, "rdma target (generic target).  waiting on wc_cq_hdl(%llu).", (uint64_t)transport_global_data.ep_cq_hdl);
                }
            } else {
                if (config.use_rdma_events) {
                    if (wr->op_state==RDMA_TARGET_INIT) {
                        log_debug(nnti_cq_debug_level, "rdma target (generic target).  waiting on mem_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
                    } else if (wr->op_state==RDMA_TARGET_NEED_ACK) {
                        log_debug(nnti_cq_debug_level, "rdma target (generic target).  waiting on wc_mem_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
                    }
                } else {
                    log_debug(nnti_cq_debug_level, "rdma target (generic target).  waiting on wc_mem_cq_hdl(%llu).", (uint64_t)transport_global_data.mem_cq_hdl);
                }
            }
            break;
        case REQUEST_BUFFER:
            log_debug(nnti_event_debug_level, "request queue.  waiting on wc_mem_cq_hdl(%llu).", (uint64_t)q_hdl->wc_mem_cq_hdl);
            break;
        case UNKNOWN_BUFFER:
        default:
            log_error(nnti_debug_level, "unknown buffer type(%llu).", gni_mem_hdl->type);
            cq_hdl=(gni_cq_handle_t)-1;
            break;
    }

    log_debug(nnti_ee_debug_level, "exit");

    return(cq_hdl);
}

static void print_failed_cq(const NNTI_buffer_t *reg_buf)
{
    gni_memory_handle        *gni_mem_hdl=NULL;
    gni_request_queue_handle *q_hdl=NULL;

    assert(reg_buf);

    q_hdl      =&transport_global_data.req_queue;
    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;

    assert(gni_mem_hdl);

    return;
}


static const NNTI_buffer_t *decode_event_buffer(
        const NNTI_buffer_t *wait_buf,
        gni_cq_entry_t      *ev_data)
{
    const NNTI_buffer_t *event_buf=NULL;
    gni_work_request  *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    if ((wait_buf != NULL) && (wait_buf->transport_private != NULL) && (((gni_memory_handle *)wait_buf->transport_private)->type == REQUEST_BUFFER)) {
        event_buf=wait_buf;

        log_debug(nnti_debug_level, "the wait buffer is a REQUEST BUFFER, so ev_data.inst_id is the request buffer index.");
    } else {
        event_buf = get_buf_bufhash((uint32_t)gni_cq_get_inst_id(*ev_data));
        if (event_buf == NULL) {
            wr = get_wr_wrhash((uint32_t)gni_cq_get_inst_id(*ev_data));
            assert(wr);
            event_buf=wr->reg_buf;
        }
        assert(event_buf);

        if (event_buf == wait_buf) {
            log_debug(nnti_debug_level, "the event buffer matches the wait buffer (ev_data.inst_id=%llu, waitbuf.buf.hash=%llu, wait_buf=%p, wr=%p, wr.hash==%llu)",
                    (uint64_t)gni_cq_get_inst_id(*ev_data), (uint64_t)hash6432shift(wait_buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf), wait_buf,
                    wr, (uint64_t)hash6432shift((uint64_t)wr));
        } else {
            log_debug(nnti_debug_level, "the event buffer does NOT match the wait buffer (ev_data.inst_id=%llu, waitbuf.buf.hash=%llu, wait_buf=%p, wr=%p, wr.hash==%llu)",
                    (uint64_t)gni_cq_get_inst_id(*ev_data), (uint64_t)hash6432shift(wait_buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf), wait_buf,
                    wr, (uint64_t)hash6432shift((uint64_t)wr));
//            log_debug(nnti_debug_level, "the ev_data.inst_id does NOT match the wait buffer (ev_data.inst_id=%lu, waitbuf.buf.hash=%lu, wait_buf=%p)",
//                    (uint64_t)gni_cq_get_inst_id(*ev_data), (uint64_t)hash6432shift(wait_buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf), wait_buf);
        }
    }

    log_debug(nnti_debug_level, "exit (event_buf==%p)", event_buf);

    return(event_buf);
}

static int process_event(
        const NNTI_buffer_t      *reg_buf,
        gni_cq_handle_t           cq_hdl,
        gni_cq_entry_t           *ev_data)
{
    int rc=NNTI_OK;
    gni_memory_handle        *gni_mem_hdl=NULL;
//    nnti_gni_work_completion *wc=NULL;
    gni_work_request         *wr=NULL;

    const NNTI_buffer_t      *event_buf=NULL;

    log_level debug_level=nnti_debug_level;

    log_debug(nnti_ee_debug_level, "enter (reg_buf=%p)", reg_buf);

    if (!GNI_CQ_STATUS_OK(*ev_data)) {
        return NNTI_EIO;
    }

    log_debug(debug_level, "ev_data.inst_id==%llu, reg_buf.buf.hash==%llu",
            (uint64_t)gni_cq_get_inst_id(*ev_data),
            (uint64_t)hash6432shift(reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf));

    event_buf=reg_buf;
    gni_mem_hdl=(gni_memory_handle *)event_buf->transport_private;
    assert(gni_mem_hdl);

    wr=first_incomplete_wr(gni_mem_hdl);
    assert(wr);

    log_debug(debug_level, "event_buf=%p; wr=%p; wr->last_op=%lld; wr.hash==%llu", event_buf, wr, (int64_t)wr->last_op, (uint64_t)hash6432shift((uint64_t)wr));
    log_debug(debug_level, "gni_mem_hdl->type==%llu", (uint64_t)gni_mem_hdl->type);

//    print_wc(&wr->wc);

    debug_level=nnti_debug_level;
    switch (gni_mem_hdl->type) {
        case SEND_BUFFER:
            if (wr->last_op==GNI_OP_SEND_REQUEST) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(debug_level, "SEND request event - event_buf==%p, op_state==%d", event_buf, wr->op_state);

                    log_debug(debug_level, "calling GetComplete(fma put send request)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(fma send request post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

                    log_debug(debug_level, "SEND request completion - event_buf==%p", event_buf);
                    wr->op_state=RDMA_WRITE_NEED_ACK;
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    log_debug(debug_level, "SEND request Work Completion completion - event_buf==%p", event_buf);

                    log_debug(debug_level, "calling GetComplete(send request)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(fma send request post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

//                    memcpy(wc, &gni_mem_hdl->wc, sizeof(nnti_gni_work_completion));

                    wr->op_state=SEND_COMPLETE;
//                    wc->byte_len   =wr->wc.byte_len;
                    wr->wc.byte_offset=wr->wc.src_offset;
                }
            } else if (wr->last_op==GNI_OP_SEND_BUFFER) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(debug_level, "SEND buffer event - event_buf==%p, op_state==%d", event_buf, wr->op_state);

                    log_debug(debug_level, "calling GetComplete(fma put send buffer)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(fma send buffer post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

                    log_debug(debug_level, "SEND buffer completion - event_buf==%p", event_buf);
                    wr->op_state=RDMA_WRITE_NEED_ACK;
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    log_debug(debug_level, "SEND buffer Work Completion completion - event_buf==%p", event_buf);

                    log_debug(debug_level, "calling GetComplete(send buffer ACK)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(fma send buffer ACK post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

                    wr->op_state = RDMA_WRITE_COMPLETE;
//                    wc->byte_len   =wr->wc.byte_len;
                    wr->wc.byte_offset=wr->wc.src_offset;
                }
            } else if (wr->last_op==GNI_OP_PUT_INITIATOR) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(debug_level, "RDMA write event - event_buf==%p, op_state==%d", event_buf, wr->op_state);

                    log_debug(debug_level, "calling GetComplete(fma put send)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(fma put send post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

                    log_debug(debug_level, "RDMA write (initiator) completion - event_buf==%p", event_buf);
                    if (config.use_rdma_target_ack) {
                        wr->op_state=RDMA_WRITE_NEED_ACK;
                    } else {
                        wr->op_state = RDMA_WRITE_COMPLETE;
                        //                    wc->byte_len   =wr->wc.byte_len;
                        wr->wc.byte_offset=wr->wc.src_offset;
                    }
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    log_debug(debug_level, "RDMA write ACK (initiator) completion - event_buf==%p", event_buf);

                    log_debug(debug_level, "calling GetComplete(fma put send ACK)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(fma put send ACK post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

                    wr->op_state = RDMA_WRITE_COMPLETE;
//                    wc->byte_len   =wr->wc.byte_len;
                    wr->wc.byte_offset=wr->wc.src_offset;
                }
            }
            break;
        case PUT_SRC_BUFFER:
            wr->last_op=GNI_OP_PUT_INITIATOR;
            if (config.use_rdma_events) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(debug_level, "RDMA write event - event_buf==%p, op_state==%d", event_buf, wr->op_state);

                    log_debug(debug_level, "calling GetComplete(rdma put src)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(put src post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

                    log_debug(debug_level, "RDMA write (initiator) completion - event_buf==%p", event_buf);
                    if (config.use_rdma_target_ack) {
                        wr->op_state=RDMA_WRITE_NEED_ACK;
                    } else {
                        wr->op_state = RDMA_WRITE_COMPLETE;
                        //                wc->byte_len   =wr->wc.byte_len;
                        wr->wc.byte_offset=wr->wc.src_offset;
                    }
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    log_debug(debug_level, "RDMA write ACK (initiator) completion - event_buf==%p", event_buf);

                    log_debug(debug_level, "calling GetComplete(rdma put src ACK)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(put src ACK post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

                    wr->op_state = RDMA_WRITE_COMPLETE;
                    //                wc->byte_len   =wr->wc.byte_len;
                    wr->wc.byte_offset=wr->wc.src_offset;
                }
            } else {
                log_debug(debug_level, "RDMA write ACK (initiator) completion - event_buf==%p", event_buf);

                log_debug(debug_level, "calling GetComplete(rdma put src ACK)");
                rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(put src ACK post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                print_post_desc(wr->post_desc_ptr);

                wr->op_state = RDMA_WRITE_COMPLETE;
                //            wc->byte_len   =wr->wc.byte_len;
                wr->wc.byte_offset=wr->wc.src_offset;
            }
            break;
        case GET_DST_BUFFER:
            wr->last_op=GNI_OP_GET_INITIATOR;
            if (config.use_rdma_events) {
                if (wr->op_state==RDMA_READ_INIT) {
                    log_debug(debug_level, "RDMA read event - event_buf==%p, op_state==%d", event_buf, wr->op_state);

                    log_debug(debug_level, "calling GetComplete(rdma get dst)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(get dst post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

                    log_debug(debug_level, "RDMA read (initiator) completion - event_buf==%p", event_buf);
                    if (config.use_rdma_target_ack) {
                        wr->op_state=RDMA_READ_NEED_ACK;
                    } else {
                        wr->op_state = RDMA_READ_COMPLETE;
                        //                wc->byte_len   =wr->wc.byte_len;
                        wr->wc.byte_offset=wr->wc.dest_offset;
                    }
                } else if (wr->op_state==RDMA_READ_NEED_ACK) {
                    log_debug(debug_level, "RDMA read ACK (initiator) completion - event_buf==%p", event_buf);

                    log_debug(debug_level, "calling GetComplete(rdma get dst ACK)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(get dst ACK post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

                    wr->op_state = RDMA_READ_COMPLETE;
                    //                wc->byte_len   =wr->wc.byte_len;
                    wr->wc.byte_offset=wr->wc.dest_offset;
                }
            } else {
                log_debug(debug_level, "RDMA read ACK (initiator) completion - event_buf==%p", event_buf);

                log_debug(debug_level, "calling GetComplete(rdma get dst ACK)");
                rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(get dst ACK post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                print_post_desc(wr->post_desc_ptr);

                wr->op_state = RDMA_READ_COMPLETE;
                //            wc->byte_len   =wr->wc.byte_len;
                wr->wc.byte_offset=wr->wc.dest_offset;
            }
            break;
        case REQUEST_BUFFER:
            {
            uint64_t  index=0;
            gni_request_queue_handle *q=&transport_global_data.req_queue;

            wr->last_op=GNI_OP_NEW_REQUEST;

            if (q->wc_buffer[q->req_processed].ack_received==1) {
                /* requests came out of order.  cleanup out of order requests here. */

                nnti_gni_work_completion *tmp_wc=&q->wc_buffer[q->req_processed];
                GNI_CQ_SET_INST_ID(*ev_data, tmp_wc->inst_id);

                log_debug(debug_level, "recv completion - event_buf=%p processing=%llu", event_buf, q->req_processed);

                wr->wc=q->wc_buffer[q->req_processed];

                wr->op_state = RECV_COMPLETE;

                q->wc_buffer[q->req_processed].ack_received=0;

                q->req_processed++;
                q->total_req_processed++;

            } else {
                index = (uint64_t)gni_cq_get_inst_id(*ev_data);
                log_debug(debug_level, "wc_index(%llu)", index);
                nnti_gni_work_completion *tmp_wc=&q->wc_buffer[index];
                tmp_wc->ack_received=1;
                GNI_CQ_SET_INST_ID(*ev_data, tmp_wc->inst_id);

                if ((q->req_processed < q->req_count) &&
                    (q->wc_buffer[q->req_processed].ack_received==0)) {

                    log_debug(nnti_event_debug_level, "request received out of order (received index(%llu) ; received ack(%llu) ; waiting(%llu)) ; waiting ack(%llu)",
                            index, q->wc_buffer[index].ack_received, q->req_processed, q->wc_buffer[q->req_processed].ack_received);
                } else {
                    log_debug(debug_level, "recv completion - event_buf=%p processing=%llu", event_buf, q->req_processed);

                    wr->wc=q->wc_buffer[q->req_processed];

                    wr->op_state = RECV_COMPLETE;

                    q->wc_buffer[q->req_processed].ack_received=0;

                    q->req_processed++;
                    q->total_req_processed++;
                }
            }

            if (q->req_processed > q->req_count) {
                log_error(debug_level, "req_processed(%llu) > req_count(%llu) fail",
                        q->req_processed, q->req_count);
            }
            if (q->req_processed == (q->req_count/2)) {
                if (q->req_index >= q->req_processed) {
                    reset_req_index(q);
                    send_unblock(q);
                    log_debug(nnti_event_debug_level, "resetting req_processed(%llu) total_req_processed(%llu)",
                            q->req_processed, q->total_req_processed);
                } else {
                    log_debug(nnti_event_debug_level, "skipping reset req_processed(%llu) total_req_processed(%llu)",
                            q->req_processed, q->total_req_processed);
                }
            }
            if (q->req_processed == q->req_processed_reset_limit) {
                q->req_processed=0;
            }
            log_debug(nnti_event_debug_level, "current req_processed(%llu) req_count(%llu)",
                    q->req_processed, q->req_count);
            log_debug(nnti_event_debug_level, "current received index(%llu) ; received ack(%llu) ; waiting(%llu)) ; waiting ack(%llu)",
                    index, q->wc_buffer[index].ack_received, q->req_processed, q->wc_buffer[q->req_processed].ack_received);
            }
            break;
        case RECEIVE_BUFFER:
            wr->last_op=GNI_OP_PUT_TARGET;
            if (config.use_rdma_events) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(debug_level, "RDMA write event - event_buf==%p, op_state==%d", event_buf, wr->op_state);
                    log_debug(debug_level, "RDMA write (receive buffer) completion - event_buf==%p", event_buf);
                    wr->op_state=RDMA_WRITE_NEED_ACK;
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    log_debug(debug_level, "RDMA write ACK (receive buffer) completion - event_buf==%p", event_buf);

                    wr->op_state = RDMA_WRITE_COMPLETE;
                    //                wc->byte_len   =wr->wc.byte_len;
                    wr->wc.byte_offset=wr->wc.dest_offset;
                }
            } else {
                log_debug(debug_level, "RDMA write ACK (receive buffer) completion - event_buf==%p", event_buf);

                wr->op_state = RDMA_WRITE_COMPLETE;
                //            wc->byte_len   =wr->wc.byte_len;
                wr->wc.byte_offset=wr->wc.dest_offset;
            }
            break;
        case PUT_DST_BUFFER:
            wr->last_op=GNI_OP_PUT_TARGET;
            if (config.use_rdma_events) {
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(debug_level, "RDMA write event - event_buf==%p, op_state==%d", event_buf, wr->op_state);
                    log_debug(debug_level, "RDMA write (target) completion - event_buf==%p", event_buf);
                    if (config.use_rdma_target_ack) {
                        wr->op_state=RDMA_WRITE_NEED_ACK;
                    } else {
                        wr->op_state = RDMA_WRITE_COMPLETE;
                        //                wc->byte_len   =wr->wc.byte_len;
                        wr->wc.byte_offset=wr->wc.dest_offset;
                    }
                } else if (wr->op_state==RDMA_WRITE_NEED_ACK) {
                    log_debug(debug_level, "RDMA write ACK (target) completion - event_buf==%p", event_buf);

                    wr->op_state = RDMA_WRITE_COMPLETE;
                    //                wc->byte_len   =wr->wc.byte_len;
                    wr->wc.byte_offset=wr->wc.dest_offset;
                }
            } else {
                log_debug(debug_level, "RDMA write ACK (target) completion - event_buf==%p", event_buf);

                wr->op_state = RDMA_WRITE_COMPLETE;
                //            wc->byte_len   =wr->wc.byte_len;
                wr->wc.byte_offset=wr->wc.dest_offset;
            }
            break;
        case GET_SRC_BUFFER:
            wr->last_op=GNI_OP_GET_TARGET;
//            if (config.rdma_mode==RDMA_CROSSOVER) {
//                if (wr->op_state==RDMA_READ_INIT) {
//                    log_debug(debug_level, "RDMA read event - event_buf==%p, op_state==%d", event_buf, wr->op_state);
//                    log_debug(debug_level, "RDMA read ACK (target) completion - event_buf==%p", event_buf);
//
//                    wr->op_state = RDMA_READ_COMPLETE;
//                    //                wc->byte_len   =wr->wc.byte_len;
//                    wr->wc.byte_offset=wr->wc.src_offset;
//                }
//            } else {
                if (config.use_rdma_events) {
                    if (wr->op_state==RDMA_READ_INIT) {
                        log_debug(debug_level, "RDMA read event - event_buf==%p, op_state==%d", event_buf, wr->op_state);
                        log_debug(debug_level, "RDMA read (target) completion - event_buf==%p", event_buf);
                        if (config.use_rdma_target_ack) {
                            wr->op_state=RDMA_READ_NEED_ACK;
                        } else {
                            wr->op_state = RDMA_READ_COMPLETE;
                            //                wc->byte_len   =wr->wc.byte_len;
                            wr->wc.byte_offset=wr->wc.src_offset;
                        }
                    } else if (wr->op_state==RDMA_READ_NEED_ACK) {
                        log_debug(debug_level, "RDMA read ACK (target) completion - event_buf==%p", event_buf);

                        wr->op_state = RDMA_READ_COMPLETE;
                        //                wc->byte_len   =wr->wc.byte_len;
                        wr->wc.byte_offset=wr->wc.src_offset;
                    }
                } else {
                    log_debug(debug_level, "RDMA read ACK (target) completion - event_buf==%p", event_buf);

                    wr->op_state = RDMA_READ_COMPLETE;
                    //            wc->byte_len   =wr->wc.byte_len;
                    wr->wc.byte_offset=wr->wc.src_offset;
                }
//            }
            break;
        case RDMA_TARGET_BUFFER:
            log_debug(debug_level, "RDMA target completion - event_buf=%p, op_state=%d, last_op=%d",
                    event_buf, wr->op_state, wr->last_op);
            if ((wr->last_op==GNI_OP_GET_INITIATOR) ||
                (wr->last_op==GNI_OP_PUT_INITIATOR)) {

                if (wr->op_state==RDMA_TARGET_INIT) {
                    log_debug(debug_level, "RDMA target event - event_buf==%p, op_state==%d", event_buf, wr->op_state);

                    log_debug(debug_level, "calling GetComplete(rdma target (mem cq)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(post_desc_ptr) failed: %d", rc);
                    print_post_desc(wr->post_desc_ptr);

                    if (config.use_rdma_target_ack) {
                        wr->op_state=RDMA_TARGET_NEED_ACK;
                    } else {
                        wr->op_state = RDMA_TARGET_COMPLETE;
                        //                    wc->byte_len   =wr->wc.byte_len;
                        wr->wc.byte_offset=wr->wc.dest_offset;
                    }
                } else if (wr->op_state==RDMA_TARGET_NEED_ACK) {

                    log_debug(debug_level, "RDMA target ACK completion - event_buf==%p", event_buf);

                    log_debug(debug_level, "calling GetComplete(rdma target ACK)");
                    rc=GNI_GetCompleted (cq_hdl, *ev_data, &wr->post_desc_ptr);
                    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "GetCompleted(rdma target ACK post_desc_ptr(%p)) failed: %d", wr->post_desc_ptr, rc);
                    print_post_desc(wr->post_desc_ptr);

                    wr->op_state = RDMA_TARGET_COMPLETE;
//                    wc->byte_len   =wr->wc.byte_len;
                    wr->wc.byte_offset=wr->wc.dest_offset;
                }
            } else {
//                if (config.rdma_mode==RDMA_CROSSOVER) {
//                    if (wr->op_state==RDMA_TARGET_INIT) {
//
//                        log_debug(debug_level, "RDMA target event - event_buf==%p, op_state==%d", event_buf, wr->op_state);
//                        log_debug(debug_level, "RDMA target completion - event_buf==%p", event_buf);
//
//                        wr->op_state = RDMA_TARGET_COMPLETE;
//                        //                    wc->byte_len   =wr->wc.byte_len;
//                        wr->wc.byte_offset=wr->wc.dest_offset;
//                    }
//                } else {
                    if (config.use_rdma_events) {
                        if (wr->op_state==RDMA_TARGET_INIT) {

                            log_debug(debug_level, "RDMA target event - event_buf==%p, op_state==%d", event_buf, wr->op_state);
                            if (config.use_rdma_target_ack) {
                                wr->op_state=RDMA_TARGET_NEED_ACK;
                            } else {
                                wr->op_state = RDMA_TARGET_COMPLETE;
                                //                    wc->byte_len   =wr->wc.byte_len;
                                wr->wc.byte_offset=wr->wc.dest_offset;
                            }
                        } else if (wr->op_state==RDMA_TARGET_NEED_ACK) {

                            log_debug(debug_level, "RDMA target completion - event_buf==%p", event_buf);

                            wr->op_state = RDMA_TARGET_COMPLETE;
                            //                    wc->byte_len   =wr->wc.byte_len;
                            wr->wc.byte_offset=wr->wc.dest_offset;
                        }
                    } else {
                        log_debug(debug_level, "RDMA target completion - event_buf==%p", event_buf);

                        wr->op_state = RDMA_TARGET_COMPLETE;
                        //                wc->byte_len   =wr->wc.byte_len;
                        wr->wc.byte_offset=wr->wc.dest_offset;
                    }
//                }
            }
            break;
    }

    print_wc(&wr->wc);

//    gni_mem_hdl->last_wc   = *wc;
//
//    switch (gni_mem_hdl->type) {
//        case RECEIVE_BUFFER:
//        case PUT_DST_BUFFER:
//        case GET_SRC_BUFFER:
//        case RDMA_TARGET_BUFFER:
//            gni_mem_hdl->last_peer = wc->inst_id;
//            break;
//    }

    log_debug(nnti_ee_debug_level, "exit");
    return (rc);
}

static NNTI_result_t post_recv_work_request(
        NNTI_buffer_t *reg_buf)
{
    gni_work_request *wr=NULL;
    gni_memory_handle *gni_mem_hdl=NULL;

    log_debug(nnti_ee_debug_level, "enter (reg_buf=%p)", reg_buf);

    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    if (config.use_wr_pool) {
        if (need_wc_mem_cq(gni_mem_hdl) == 1) {
            wr=wr_pool_target_pop();
        } else {
            wr=wr_pool_initiator_pop();
        }
    } else {
        wr=(gni_work_request *)calloc(1, sizeof(gni_work_request));
    }
    assert(wr);

    wr->reg_buf = reg_buf;

    if (!config.use_wr_pool) {
        register_wc(wr);
    }

    if (gni_mem_hdl->type==RECEIVE_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (gni_mem_hdl->type==GET_SRC_BUFFER) {
        wr->op_state=RDMA_READ_INIT;

    } else if (gni_mem_hdl->type==PUT_DST_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (gni_mem_hdl->type==RDMA_TARGET_BUFFER) {
        wr->op_state=RDMA_TARGET_INIT;

    }

    nthread_lock(&gni_mem_hdl->wr_queue_lock);
    gni_mem_hdl->wr_queue.push_back(wr);
    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

    reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.wc_addr           = (uint64_t)&wr->wc;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword1 = wr->wc_mem_hdl.qword1;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword2 = wr->wc_mem_hdl.qword2;

    log_debug(nnti_ee_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t repost_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        gni_work_request *wr)
{
    gni_memory_handle *gni_mem_hdl=NULL;

    log_debug(nnti_ee_debug_level, "enter (reg_buf=%p)", reg_buf);

    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    memset(&wr->wc, 0, sizeof(nnti_gni_work_completion));
    wr->last_op       =0;
    wr->is_op_complete=FALSE;
    wr->op_state      =BUFFER_INIT;

    if (gni_mem_hdl->type==RECEIVE_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (gni_mem_hdl->type==GET_SRC_BUFFER) {
        wr->op_state=RDMA_READ_INIT;

    } else if (gni_mem_hdl->type==PUT_DST_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (gni_mem_hdl->type==RDMA_TARGET_BUFFER) {
        wr->op_state=RDMA_TARGET_INIT;

    }

    gni_mem_hdl->wr_queue.push_back(wr);

    reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.wc_addr           = (uint64_t)&wr->wc;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword1 = wr->wc_mem_hdl.qword1;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword2 = wr->wc_mem_hdl.qword2;

    log_debug(nnti_ee_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static int8_t is_wr_complete(
        const gni_work_request *wr)
{
    int8_t rc=FALSE;
    gni_memory_handle *gni_mem_hdl=NULL;

    gni_mem_hdl=(gni_memory_handle *)wr->reg_buf->transport_private;
    assert(gni_mem_hdl);

    switch (gni_mem_hdl->type) {
        case SEND_BUFFER:
            if ((wr->op_state == SEND_COMPLETE) ||
                (wr->op_state == RDMA_WRITE_COMPLETE)) {
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
            if (wr->op_state == RECV_COMPLETE) {
                rc=TRUE;
            }
            break;
        case RECEIVE_BUFFER:
            if (wr->op_state == RDMA_WRITE_COMPLETE) {
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

    log_debug(nnti_ee_debug_level, "exit (rc=%d)", rc);
    return(rc);
}

static gni_work_request *first_incomplete_wr(
        gni_memory_handle *gni_mem_hdl)
{
    gni_work_request  *wr=NULL;

    log_debug(nnti_ee_debug_level, "enter");

    assert(gni_mem_hdl);

    nthread_lock(&gni_mem_hdl->wr_queue_lock);

    if (gni_mem_hdl->wr_queue.empty()) {
        log_debug(nnti_ee_debug_level, "work request queue is empty");
    } else {
        wr_queue_iter_t i;
        for (i=gni_mem_hdl->wr_queue.begin(); i != gni_mem_hdl->wr_queue.end(); i++) {
            wr=*i;
            assert(wr);
            if (is_wr_complete(wr) == FALSE) {
                log_debug(nnti_ee_debug_level, "incomplete work request found");
                break;
            }
        }
        if (wr==NULL) {
            log_debug(nnti_ee_debug_level, "all work requests complete");
        }
    }

    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

    log_debug(nnti_ee_debug_level, "exit (wr=%p)", wr);
    return(wr);
}

static int8_t is_wr_queue_empty(
        const NNTI_buffer_t *reg_buf)
{
    int8_t rc=FALSE;
    gni_memory_handle *gni_mem_hdl=NULL;

    log_debug(nnti_ee_debug_level, "enter");

    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    if (gni_mem_hdl->wr_queue.empty()) {
        rc=TRUE;
    }

    log_debug(nnti_ee_debug_level, "exit (rc=%d)", rc);
    return(rc);
}

static int8_t is_buf_op_complete(
        const NNTI_buffer_t *reg_buf)
{
    int8_t rc=FALSE;
    gni_memory_handle *gni_mem_hdl=NULL;
    gni_work_request  *wr=NULL;

    log_debug(nnti_ee_debug_level, "enter");

    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    nthread_lock(&gni_mem_hdl->wr_queue_lock);

    if (is_wr_queue_empty(reg_buf) == TRUE) {
        log_debug(nnti_ee_debug_level, "work request queue is empty - return FALSE");
        rc=FALSE;
    } else {
        wr=gni_mem_hdl->wr_queue.front();
        assert(wr);

        rc = is_wr_complete(wr);
    }

    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

    log_debug(nnti_ee_debug_level, "exit (rc=%d)", rc);
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
        int                   nnti_rc,
        gni_cq_entry_t       *ev_data,
        NNTI_status_t        *status)
{
    gni_connection    *conn       =NULL;
    gni_memory_handle *gni_mem_hdl=NULL;
    gni_work_request  *wr=NULL;

    memset(status, 0, sizeof(NNTI_status_t));
    status->op     = remote_op;
    status->result = (NNTI_result_t)nnti_rc;
    if (nnti_rc==NNTI_OK) {
        gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
        assert(gni_mem_hdl);
        wr=gni_mem_hdl->wr_queue.front();
        assert(wr);

        if (gni_mem_hdl->type == REQUEST_BUFFER) {
            conn = get_conn_instance(GNI_CQ_GET_INST_ID(*ev_data));
        } else {
            if (wr->peer_instance) {
                conn = get_conn_instance(wr->peer_instance);
            } else {
                conn = get_conn_instance(wr->wc.inst_id);
            }
        }
        assert(conn);

        status->start  = (uint64_t)reg_buf->payload;
        status->offset = wr->wc.byte_offset;
        status->length = wr->wc.byte_len;
        switch (wr->last_op) {
            case GNI_OP_PUT_INITIATOR:
            case GNI_OP_GET_TARGET:
            case GNI_OP_SEND_REQUEST:
            case GNI_OP_SEND_BUFFER:
                create_peer(&status->src,
                        transport_global_data.listen_name,
                        transport_global_data.listen_addr,
                        transport_global_data.listen_port,
                        transport_global_data.alps_info.ptag,
                        transport_global_data.alps_info.cookie,
                        transport_global_data.instance);
                create_peer(&status->dest,
                        conn->peer_name,
                        conn->peer_addr,
                        conn->peer_port,
                        conn->peer_ptag,
                        conn->peer_cookie,
                        conn->peer_instance);
                break;
            case GNI_OP_GET_INITIATOR:
            case GNI_OP_PUT_TARGET:
            case GNI_OP_NEW_REQUEST:
            case GNI_OP_RECEIVE:
                create_peer(&status->src,
                        conn->peer_name,
                        conn->peer_addr,
                        conn->peer_port,
                        conn->peer_ptag,
                        conn->peer_cookie,
                        conn->peer_instance);
                create_peer(&status->dest,
                        transport_global_data.listen_name,
                        transport_global_data.listen_addr,
                        transport_global_data.listen_port,
                        transport_global_data.alps_info.ptag,
                        transport_global_data.alps_info.cookie,
                        transport_global_data.instance);
                break;
            }
    }
}

static void create_peer(NNTI_peer_t *peer, char *name, NNTI_ip_addr addr, NNTI_tcp_port port, uint32_t ptag, uint32_t cookie, NNTI_instance_id instance)
{
    log_debug(nnti_ee_debug_level, "enter");

    sprintf(peer->url, "gni://%s:%u/?ptag=%llu&cookie=%llu", name, ntohs(port), (uint64_t)ptag, (uint64_t)cookie);

    peer->peer.transport_id                       =NNTI_TRANSPORT_GEMINI;
    peer->peer.NNTI_remote_process_t_u.gni.addr   =addr;
    peer->peer.NNTI_remote_process_t_u.gni.port   =port;
    peer->peer.NNTI_remote_process_t_u.gni.inst_id=instance;

    log_debug(nnti_ee_debug_level, "exit");
}

static void copy_peer(NNTI_peer_t *src, NNTI_peer_t *dest)
{
    log_debug(nnti_ee_debug_level, "enter");

    strncpy(dest->url, src->url, NNTI_URL_LEN);
    dest->url[NNTI_URL_LEN-1]='\0';

    dest->peer.transport_id                       =NNTI_TRANSPORT_GEMINI;
    dest->peer.NNTI_remote_process_t_u.gni.addr   =src->peer.NNTI_remote_process_t_u.gni.addr;
    dest->peer.NNTI_remote_process_t_u.gni.port   =src->peer.NNTI_remote_process_t_u.gni.port;
    dest->peer.NNTI_remote_process_t_u.gni.inst_id=src->peer.NNTI_remote_process_t_u.gni.inst_id;

    log_debug(nnti_ee_debug_level, "exit");
}

static int init_server_listen_socket()
{
    NNTI_result_t rc=NNTI_OK;
    trios_declare_timer(call_time);
    int flags;
    struct hostent *host_entry;
    struct sockaddr_in skin;
    socklen_t skin_size=sizeof(struct sockaddr_in);

    trios_start_timer(call_time);
    transport_global_data.listen_sock = socket(AF_INET, SOCK_STREAM, 0);
    trios_stop_timer("socket", call_time);
    if (transport_global_data.listen_sock < 0)
        log_error(nnti_debug_level, "failed to create tcp socket: %s", strerror(errno));

    flags = 1;
    trios_start_timer(call_time);
    if (setsockopt(transport_global_data.listen_sock, SOL_SOCKET, SO_REUSEADDR, &flags, sizeof(flags)) < 0)
        log_error(nnti_debug_level, "failed to set tcp socket REUSEADDR flag: %s", strerror(errno));
    trios_stop_timer("setsockopt", call_time);

    flags=fcntl(transport_global_data.listen_sock, F_GETFL, 0);
    fcntl(transport_global_data.listen_sock, F_SETFL, flags | O_NONBLOCK);

    log_debug(nnti_debug_level, "listen_name (%s).", transport_global_data.listen_name);
    if (transport_global_data.listen_name[0]!='\0') {
        log_debug(nnti_debug_level, "using hostname from command-line (%s).", transport_global_data.listen_name);
    } else {
        trios_start_timer(call_time);
        gethostname(transport_global_data.listen_name, NNTI_HOSTNAME_LEN);
        trios_stop_timer("gethostname", call_time);
        log_debug(nnti_debug_level, "hostname not given on command-line.  using gethostname() result (%s).", transport_global_data.listen_name);
    }

    /* lookup the host provided on the command line */
    trios_start_timer(call_time);
    host_entry = gethostbyname(transport_global_data.listen_name);
    trios_stop_timer("gethostbyname", call_time);
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
    trios_start_timer(call_time);
    if (bind(transport_global_data.listen_sock, (struct sockaddr *)&skin, skin_size) < 0) {
        if (errno == EINTR) {
            goto retry;
        } else {
            log_error(nnti_debug_level, "failed to bind tcp socket: %s", strerror(errno));
        }
    }
    trios_stop_timer("bind", call_time);


    /* after the bind, get the "name" for the socket.  the "name" contains the port assigned by the kernel. */
    trios_start_timer(call_time);
    getsockname(transport_global_data.listen_sock, (struct sockaddr *)&skin, &skin_size);
    trios_stop_timer("getsockname", call_time);
    transport_global_data.listen_addr = (uint32_t)skin.sin_addr.s_addr;
    transport_global_data.listen_port = (uint16_t)skin.sin_port;
    log_debug(nnti_debug_level, "listening on ip(%s) addr(%u) port(%u)",
            transport_global_data.listen_name,
            (unsigned int)ntohl(skin.sin_addr.s_addr),
            (unsigned int)ntohs(skin.sin_port));
    trios_start_timer(call_time);
    if (listen(transport_global_data.listen_sock, 1024) < 0)
        log_error(nnti_debug_level, "failed to listen on tcp socket: %s", strerror(errno));
    trios_stop_timer("listen", call_time);

    return rc;
}

static void transition_connection_to_ready(
        int sock,
        gni_connection *conn)
{
    int rc=NNTI_OK;
    trios_declare_timer(callTime);

    trios_start_timer(callTime);
    /* final sychronization to ensure both sides have posted RTRs */
    rc = tcp_exchange(sock, 0, &rc, &rc, sizeof(rc));
    trios_stop_timer("transition tcp_exchange", callTime);
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
    int rc=0;

    if (is_server) {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_read(sock, incoming, len);
        trios_stop_timer("tcp_read", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "server failed to read GNI connection info: errno=%d", errno);
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
            log_warn(nnti_debug_level, "client failed to write GNI connection info: errno=%d", errno);
            goto out;
        }
    }

    if (is_server) {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_write(sock, outgoing, len);
        trios_stop_timer("tcp_write", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "server failed to write GNI connection info: errno=%d", errno);
            goto out;
        }
    } else {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_read(sock, incoming, len);
        trios_stop_timer("tcp_read", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "client failed to read GNI connection info: errno=%d", errno);
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
        gni_connection *c,
        int sock)
{
    int rc;

    /*
     * Values passed through TCP to permit Gemini connection
     */
    struct {
        char             name[NNTI_HOSTNAME_LEN];
        uint32_t         addr;
        uint32_t         port;
        NNTI_instance_id instance;
        alpsAppGni_t     alps_info;

        uint32_t         client_has_recv_queue;
    } instance_in, instance_out;
    struct {
        nnti_gni_recv_queue_attrs server_attrs;
    } sa_in, sa_out;
    struct {
        nnti_gni_recv_queue_flow_control_attrs client_attrs;
    } ca_out, ca_in;

    trios_declare_timer(call_time);

    c->connection_type=CLIENT_CONNECTION;

    /*
     * Prepare to exchange TCP and ALPS parameters with the server
     */
    memset(&instance_out, 0, sizeof(instance_out));
    strcpy(instance_out.name, transport_global_data.listen_name);
    instance_out.addr      = htonl((uint32_t)transport_global_data.listen_addr);
    instance_out.port      = htonl((uint32_t)transport_global_data.listen_port);
    instance_out.instance  = htonl(transport_global_data.instance);
    instance_out.alps_info = transport_global_data.alps_info;
    if ((transport_global_data.req_queue.req_size  > 0) &&
        (transport_global_data.req_queue.req_count > 0)) {
        instance_out.client_has_recv_queue = 1;
    }

    /*
     * Exchange TCP and ALPS parameters with the server
     */
    trios_start_timer(call_time);
    rc = tcp_exchange(sock, 0, &instance_in, &instance_out, sizeof(instance_in));
    trios_stop_timer("tcp_exchange", call_time);
    if (rc)
        goto out;

    /*
     * Record the server's parameters
     */
    c->peer_name      = strdup(instance_in.name);
    c->peer_addr      = ntohl(instance_in.addr);
    c->peer_port      = ntohl(instance_in.port);
    c->peer_instance  = ntohl(instance_in.instance);
    c->peer_alps_info = instance_in.alps_info;

    /*
     * Read the receive queue attributes from the server
     */
    memset(&sa_in, 0, sizeof(sa_in));
    trios_start_timer(call_time);
    rc = tcp_read(sock, &sa_in, sizeof(sa_in));
    trios_stop_timer("read server queue attrs", call_time);
    if (rc == sizeof(sa_in)) {
        rc=0;
    }
    if (rc)
        goto out;

    /*
     * Record the server's attributes
     */
    c->queue_remote_attrs.server=sa_in.server_attrs;

    /*
     * Setup flow control attributes
     */
    client_req_queue_init(c);

    /*
     * Write receive queue flow control attributes to the server
     */
    memset(&ca_out, 0, sizeof(ca_out));
    ca_out.client_attrs.unblock_buffer_addr=c->queue_local_attrs.unblock_buffer_addr;
    ca_out.client_attrs.unblock_mem_hdl    =c->queue_local_attrs.unblock_mem_hdl;

    trios_start_timer(call_time);
    rc = tcp_write(sock, &ca_out, sizeof(ca_out));
    trios_stop_timer("write client queue attrs", call_time);
    if (rc == sizeof(ca_out)) {
        rc=0;
    }
    if (rc)
        goto out;



    /*
     * If the client registered a receive queue (bidirectional connection)...
     */
    if (instance_out.client_has_recv_queue == 1) {
        memset(&sa_out, 0, sizeof(sa_out));

        gni_memory_handle *gni_mem_hdl=(gni_memory_handle *)transport_global_data.req_queue.reg_buf->transport_private;

        /*
         * Then write the receive queue attributes to the server
         */
        sa_out.server_attrs.req_index_addr   =transport_global_data.req_queue.req_index_addr;
        sa_out.server_attrs.req_index_mem_hdl=transport_global_data.req_queue.req_index_mem_hdl;
        sa_out.server_attrs.req_buffer_addr  =(uint64_t)transport_global_data.req_queue.req_buffer;
        sa_out.server_attrs.req_size         =transport_global_data.req_queue.req_size;
        sa_out.server_attrs.req_count        =transport_global_data.req_queue.req_count;
        sa_out.server_attrs.req_mem_hdl      =gni_mem_hdl->mem_hdl;
        sa_out.server_attrs.wc_buffer_addr   =(uint64_t)transport_global_data.req_queue.wc_buffer;
        sa_out.server_attrs.wc_mem_hdl       =transport_global_data.req_queue.wc_mem_hdl;

        trios_start_timer(call_time);
        rc = tcp_write(sock, &sa_out, sizeof(sa_out));
        trios_stop_timer("write server queue attrs", call_time);
        if (rc == sizeof(sa_out)) {
            rc=0;
        }
        if (rc)
            goto out;

        /*
         * Then read receive queue flow control attributes from the server
         */
        memset(&ca_in, 0, sizeof(ca_in));
        trios_start_timer(call_time);
        rc = tcp_read(sock, &ca_in, sizeof(ca_in));
        trios_stop_timer("read client queue attrs", call_time);
        if (rc == sizeof(ca_in)) {
            rc=0;
        }
        if (rc)
            goto out;

        /*
         * Record the server's flow control attributes
         */
        c->queue_remote_attrs.client=ca_in.client_attrs;
    }


    rc=GNI_EpCreate (c->nic_hdl, transport_global_data.ep_cq_hdl, &c->ep_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "EpCreate(c->ep_hdl) failed: %d", rc);
        goto out;
    }
    rc=GNI_EpBind (c->ep_hdl, c->peer_alps_info.local_addr, c->peer_instance);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "EpBind(c->ep_hdl) failed: %d", rc);
        goto out;
    }

out:
    return rc;
}

static int new_server_connection(
        gni_connection *c,
        int sock)
{
    int rc;

    /*
     * Values passed through TCP to permit Gemini connection
     */
    struct {
        char             name[NNTI_HOSTNAME_LEN];
        uint32_t         addr;
        uint32_t         port;
        NNTI_instance_id instance;
        alpsAppGni_t     alps_info;

        uint32_t         client_has_recv_queue;
    } instance_in, instance_out;
    struct {
        nnti_gni_recv_queue_attrs server_attrs;
    } sa_out, sa_in;
    struct {
        nnti_gni_recv_queue_flow_control_attrs client_attrs;
    } ca_in, ca_out;

    trios_declare_timer(call_time);

    assert(transport_global_data.req_queue.reg_buf);

    gni_memory_handle *gni_mem_hdl=(gni_memory_handle *)transport_global_data.req_queue.reg_buf->transport_private;


    c->connection_type=SERVER_CONNECTION;

    /*
     * Prepare to exchange TCP and ALPS parameters with the client
     */
    memset(&instance_out, 0, sizeof(instance_out));
    strcpy(instance_out.name, transport_global_data.listen_name);
    instance_out.addr      = htonl((uint32_t)transport_global_data.listen_addr);
    instance_out.port      = htonl((uint32_t)transport_global_data.listen_port);
    instance_out.instance  = htonl(transport_global_data.instance);
    instance_out.alps_info = transport_global_data.alps_info;

    /*
     * Exchange TCP and ALPS parameters with the client
     */
    trios_start_timer(call_time);
    rc = tcp_exchange(sock, 1, &instance_in, &instance_out, sizeof(instance_in));
    trios_stop_timer("tcp_exchange", call_time);
    if (rc)
        goto out;

    /*
     * Record the client's parameters
     */
    c->peer_name      = strdup(instance_in.name);
    c->peer_addr      = ntohl(instance_in.addr);
    c->peer_port      = ntohl(instance_in.port);
    c->peer_instance  = ntohl(instance_in.instance);
    c->peer_alps_info = instance_in.alps_info;
    c->peer_ptag      = instance_in.alps_info.ptag;
    c->peer_cookie    = instance_in.alps_info.cookie;
    c->cdm_hdl = transport_global_data.cdm_hdl;
    c->nic_hdl = transport_global_data.nic_hdl;

    /*
     * Write the receive queue attributes to the client
     */
    memset(&sa_out, 0, sizeof(sa_out));
    sa_out.server_attrs.req_index_addr   =transport_global_data.req_queue.req_index_addr;
    sa_out.server_attrs.req_index_mem_hdl=transport_global_data.req_queue.req_index_mem_hdl;
    sa_out.server_attrs.req_buffer_addr  =(uint64_t)transport_global_data.req_queue.req_buffer;
    sa_out.server_attrs.req_size         =transport_global_data.req_queue.req_size;
    sa_out.server_attrs.req_count        =transport_global_data.req_queue.req_count;
    sa_out.server_attrs.req_mem_hdl      =gni_mem_hdl->mem_hdl;
    sa_out.server_attrs.wc_buffer_addr   =(uint64_t)transport_global_data.req_queue.wc_buffer;
    sa_out.server_attrs.wc_mem_hdl       =transport_global_data.req_queue.wc_mem_hdl;

    trios_start_timer(call_time);
    rc = tcp_write(sock, &sa_out, sizeof(sa_out));
    trios_stop_timer("write server queue attrs", call_time);
    if (rc == sizeof(sa_out)) {
        rc=0;
    }
    if (rc)
        goto out;

    /*
     * Read receive queue flow control attributes from the client
     */
    memset(&ca_in, 0, sizeof(ca_in));
    trios_start_timer(call_time);
    rc = tcp_read(sock, &ca_in, sizeof(ca_in));
    trios_stop_timer("read client queue attrs", call_time);
    if (rc == sizeof(ca_in)) {
        rc=0;
    }
    if (rc)
        goto out;

    /*
     * Record the client's flow control attributes
     */
    c->queue_remote_attrs.client=ca_in.client_attrs;



    /*
     * If the client registered a receive queue (bidirectional connection)...
     */
    if (instance_in.client_has_recv_queue == 1) {
        /*
         * Then read the receive queue attributes from the client
         */
        memset(&sa_in, 0, sizeof(sa_in));
        trios_start_timer(call_time);
        rc = tcp_read(sock, &sa_in, sizeof(sa_in));
        trios_stop_timer("read server queue attrs", call_time);
        if (rc == sizeof(sa_in)) {
            rc=0;
        }
        if (rc)
            goto out;

        /*
         * Record the client's attributes
         */
        c->queue_remote_attrs.server=sa_in.server_attrs;

        /*
         * Setup flow control attributes
         */
        client_req_queue_init(c);

        /*
         * Then write receive queue flow control attributes to the client
         */
        memset(&ca_out, 0, sizeof(ca_out));
        ca_out.client_attrs.unblock_buffer_addr=c->queue_local_attrs.unblock_buffer_addr;
        ca_out.client_attrs.unblock_mem_hdl    =c->queue_local_attrs.unblock_mem_hdl;

        trios_start_timer(call_time);
        rc = tcp_write(sock, &ca_out, sizeof(ca_out));
        trios_stop_timer("write client queue attrs", call_time);
        if (rc == sizeof(ca_out)) {
            rc=0;
        }
        if (rc)
            goto out;
    }


    rc=GNI_EpCreate (c->nic_hdl, transport_global_data.ep_cq_hdl, &c->ep_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "EpCreate(c->ep_hdl) failed: %d", rc);
        goto out;
    }
    rc=GNI_EpBind (c->ep_hdl, c->peer_alps_info.local_addr, c->peer_instance);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "EpBind(c->ep_hdl) failed: %d", rc);
        goto out;
    }
out:
    return rc;
}

static NNTI_result_t insert_conn_peer(const NNTI_peer_t *peer, gni_connection *conn)
{
    NNTI_result_t  rc=NNTI_OK;
    addrport_key key;

    key.addr = peer->peer.NNTI_remote_process_t_u.gni.addr;
    key.port = peer->peer.NNTI_remote_process_t_u.gni.port;

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "peer",
                "insert_conn_peer", peer);
    }

    if (nthread_lock(&nnti_conn_peer_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (connections_by_peer.find(key) != connections_by_peer.end()) {
        print_peer_map();
        assert(connections_by_peer.find(key) == connections_by_peer.end());
    }
    connections_by_peer[key] = conn;   // add to connection map
    nthread_unlock(&nnti_conn_peer_lock);

    log_debug(nnti_debug_level, "peer connection added (conn=%p)", conn);

    return(rc);
}
static NNTI_result_t insert_conn_instance(const NNTI_instance_id instance, gni_connection *conn)
{
    NNTI_result_t  rc=NNTI_OK;

    if (nthread_lock(&nnti_conn_instance_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (connections_by_instance.find(instance) != connections_by_instance.end()) {
        print_instance_map();
        assert(connections_by_instance.find(instance) == connections_by_instance.end());
    }
    connections_by_instance[instance] = conn;
    nthread_unlock(&nnti_conn_instance_lock);

    log_debug(nnti_debug_level, "instance connection added (conn=%p)", conn);

    return(rc);
}
static gni_connection *get_conn_peer(const NNTI_peer_t *peer)
{
    gni_connection *conn = NULL;
    addrport_key   key;

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "looking for peer",
                "get_conn_peer", peer);
    }

    memset(&key, 0, sizeof(addrport_key));
    key.addr=peer->peer.NNTI_remote_process_t_u.gni.addr;
    key.port=peer->peer.NNTI_remote_process_t_u.gni.port;

    if (nthread_lock(&nnti_conn_peer_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (connections_by_peer.find(key) != connections_by_peer.end()) {
        conn = connections_by_peer[key];
    }
    nthread_unlock(&nnti_conn_peer_lock);

    if (conn != NULL) {
        log_debug(nnti_debug_level, "connection found");
        return conn;
    }

    log_debug(nnti_debug_level, "connection NOT found");
    print_peer_map();

    return(NULL);
}
static gni_connection *get_conn_instance(const NNTI_instance_id instance)
{
    gni_connection *conn=NULL;

    log_debug(nnti_debug_level, "looking for instance=%llu", (unsigned long long)instance);
    if (nthread_lock(&nnti_conn_instance_lock)) log_warn(nnti_debug_level, "failed to get thread lock");;
    if (connections_by_instance.find(instance) != connections_by_instance.end()) {
        conn = connections_by_instance[instance];
    }
    nthread_unlock(&nnti_conn_instance_lock);

    if (conn != NULL) {
        log_debug(nnti_debug_level, "connection found");
        return conn;
    }

    log_debug(nnti_debug_level, "connection NOT found (instance==%llu)", (uint64_t)instance);
    print_instance_map();

    return(NULL);
}
static NNTI_peer_t *get_peer_by_url(const char *url)
{
    gni_connection *conn = NULL;

    log_debug(nnti_debug_level, "looking for url=%s", url);

    conn_by_peer_iter_t i;
    if (nthread_lock(&nnti_conn_peer_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    for (i=connections_by_peer.begin(); i != connections_by_peer.end(); i++) {
        log_debug(nnti_debug_level, "peer_map key=(%llu,%llu) conn=%p, url=%s",
                (uint64_t)i->first.addr, (uint64_t)i->first.port, i->second, i->second->peer.url);
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
static gni_connection *del_conn_peer(const NNTI_peer_t *peer)
{
    gni_connection *conn=NULL;
    addrport_key    key;

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "peer",
                "get_conn_peer", peer);
    }

    memset(&key, 0, sizeof(addrport_key));
    key.addr=peer->peer.NNTI_remote_process_t_u.gni.addr;
    key.port=peer->peer.NNTI_remote_process_t_u.gni.port;

    if (nthread_lock(&nnti_conn_peer_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (connections_by_peer.find(key) != connections_by_peer.end()) {
        conn = connections_by_peer[key];
    }

    if (conn != NULL) {
        log_debug(nnti_debug_level, "connection found");
        connections_by_peer.erase(key);
        del_conn_instance(conn->peer_instance);
    } else {
        log_debug(nnti_debug_level, "connection NOT found");
    }
    nthread_unlock(&nnti_conn_peer_lock);

    return(conn);
}
static gni_connection *del_conn_instance(const NNTI_instance_id instance)
{
    gni_connection *conn=NULL;
    log_level debug_level = nnti_debug_level;

    if (nthread_lock(&nnti_conn_instance_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (connections_by_instance.find(instance) != connections_by_instance.end()) {
        conn = connections_by_instance[instance];
    }

    if (conn != NULL) {
        log_debug(debug_level, "connection found");
        connections_by_instance.erase(instance);
    } else {
        log_debug(debug_level, "connection NOT found");
    }
    nthread_unlock(&nnti_conn_instance_lock);

    return(conn);
}
static void print_peer_map()
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    conn_by_peer_iter_t i;
    for (i=connections_by_peer.begin(); i != connections_by_peer.end(); i++) {
        log_debug(nnti_debug_level, "peer_map key=(%llu,%llu) conn=%p",
                (uint64_t)i->first.addr, (uint64_t)i->first.port, i->second);
    }
}
static void print_instance_map()
{
//    if (!logging_debug(nnti_debug_level)) {
//        return;
//    }

    conn_by_inst_iter_t i;
    for (i=connections_by_instance.begin(); i != connections_by_instance.end(); i++) {
        log_debug(nnti_debug_level, "instance_map key=%llu conn=%p", i->first, i->second);
    }
}


static NNTI_result_t insert_buf_bufhash(NNTI_buffer_t *buf)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t h=hash6432shift((uint64_t)buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf);

    if (nthread_lock(&nnti_buf_bufhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");;
    assert(buffers_by_bufhash.find(h) == buffers_by_bufhash.end());
    buffers_by_bufhash[h] = buf;
    nthread_unlock(&nnti_buf_bufhash_lock);

    log_debug(nnti_debug_level, "bufhash buffer added (buf=%p)", buf);

    return(rc);
}
static NNTI_buffer_t *get_buf_bufhash(const uint32_t bufhash)
{
    NNTI_buffer_t *buf=NULL;

    log_debug(nnti_debug_level, "looking for bufhash=%llu", (uint64_t)bufhash);
    if (nthread_lock(&nnti_buf_bufhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");;
    if (buffers_by_bufhash.find(bufhash) != buffers_by_bufhash.end()) {
        buf = buffers_by_bufhash[bufhash];
    }
    nthread_unlock(&nnti_buf_bufhash_lock);

    if (buf != NULL) {
        log_debug(nnti_debug_level, "buffer found (buf=%p)", buf);
        return buf;
    }

    log_debug(nnti_debug_level, "buffer NOT found");
    print_bufhash_map();

    return(NULL);
}
static NNTI_buffer_t *del_buf_bufhash(NNTI_buffer_t *buf)
{
    uint32_t h=hash6432shift((uint64_t)buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf);
    log_level debug_level = nnti_debug_level;

    if (nthread_lock(&nnti_buf_bufhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
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

    buf_by_bufhash_iter_t i;
    for (i=buffers_by_bufhash.begin(); i != buffers_by_bufhash.end(); i++) {
        log_debug(nnti_debug_level, "bufhash_map key=%llu buf=%p", i->first, i->second);
    }
}

static NNTI_result_t insert_wr_wrhash(gni_work_request *wr)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t h=hash6432shift((uint64_t)wr);

    if (nthread_lock(&nnti_wr_wrhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    assert(wr_by_wrhash.find(h) == wr_by_wrhash.end());
    wr_by_wrhash[h] = wr;
    nthread_unlock(&nnti_wr_wrhash_lock);

    log_debug(nnti_debug_level, "wrhash work request added (wr=%p ; wr.hash=%llu)", wr, (uint64_t)h);

    return(rc);
}
static gni_work_request *get_wr_wrhash(const uint32_t wrhash)
{
    gni_work_request *wr=NULL;

    log_debug(nnti_debug_level, "looking for wrhash=%llu", (uint64_t)wrhash);
    if (nthread_lock(&nnti_wr_wrhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (wr_by_wrhash.find(wrhash) != wr_by_wrhash.end()) {
        wr = wr_by_wrhash[wrhash];
    }
    nthread_unlock(&nnti_wr_wrhash_lock);

    if (wr != NULL) {
        log_debug(nnti_debug_level, "work request found (wr=%p)", wr);
        return wr;
    }

    log_debug(nnti_debug_level, "work request NOT found");
    print_wrhash_map();

    return(NULL);
}
static gni_work_request *del_wr_wrhash(gni_work_request *wr)
{
    uint32_t h=hash6432shift((uint64_t)wr);
    log_level debug_level = nnti_debug_level;

    if (nthread_lock(&nnti_wr_wrhash_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (wr_by_wrhash.find(h) != wr_by_wrhash.end()) {
        wr = wr_by_wrhash[h];
    }

    if (wr != NULL) {
        log_debug(debug_level, "work request found");
        wr_by_wrhash.erase(h);
    } else {
        log_debug(debug_level, "work request NOT found");
    }
    nthread_unlock(&nnti_wr_wrhash_lock);

    return(wr);
}
static void print_wrhash_map()
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    wr_by_wrhash_iter_t i;
    for (i=wr_by_wrhash.begin(); i != wr_by_wrhash.end(); i++) {
        log_debug(nnti_debug_level, "wrhash_map key=%llu wr=%p", i->first, i->second);
    }
}

static NNTI_result_t wr_pool_register(
        gni_work_request *wr,
        gni_cq_handle_t   cq_hdl)
{
    int rc=GNI_RC_SUCCESS;

    log_debug(nnti_debug_level, "enter");

    rc=GNI_MemRegister (transport_global_data.nic_hdl,
            (uint64_t)&wr->wc,
            sizeof(nnti_gni_work_completion),
            cq_hdl,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &wr->wc_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "MemRegister(wc_mem_hdl) failed: rc=%d, %s", rc, strerror(errno));
        return(NNTI_EIO);
    }

    wr->wc_registered=TRUE;

    log_debug(nnti_debug_level, "exit");

    return(NNTI_OK);
}
static NNTI_result_t wr_pool_deregister(
        gni_work_request *wr)
{
    int rc=GNI_RC_SUCCESS;

    log_debug(nnti_debug_level, "enter");

    if (wr->wc_registered==FALSE) {
        log_debug(nnti_debug_level, "exit wr(%p) - not registered", wr);
        return(NNTI_OK);
    }

    rc=GNI_MemDeregister(transport_global_data.nic_hdl, &wr->wc_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) {
        log_error(nnti_debug_level, "MemDeregister(wc_mem_hdl) failed: rc=%d", rc);
        return(NNTI_EIO);
    }

    wr->wc_registered=FALSE;

    log_debug(nnti_debug_level, "exit");

    return(NNTI_OK);
}
static NNTI_result_t wr_pool_init(uint32_t pool_size)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;
    gni_work_request *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    for (i=0;i<pool_size;i++) {
        wr=(gni_work_request *)calloc(1, sizeof(gni_work_request));
        assert(wr);
        rc=wr_pool_register(wr, transport_global_data.mem_cq_hdl);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to register target work request: rc=%d", rc);
            goto cleanup;
        }
        wr->is_initiator=FALSE;
        wr_pool_target_push(wr);

        wr=(gni_work_request *)calloc(1, sizeof(gni_work_request));
        assert(wr);
        rc=wr_pool_register(wr, NULL);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to register initiator work request: rc=%d", rc);
            goto cleanup;
        }
        wr->is_initiator=TRUE;
        wr_pool_initiator_push(wr);
    }

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(rc);
}
static gni_work_request *wr_pool_target_pop(void)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;
    gni_work_request *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    if (nthread_lock(&nnti_wr_pool_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (!target_wr_pool.empty()) {
        wr=target_wr_pool.front();
        target_wr_pool.pop_front();
    }
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return(wr);
}
static gni_work_request *wr_pool_initiator_pop(void)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;
    gni_work_request *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    if (nthread_lock(&nnti_wr_pool_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    if (!initiator_wr_pool.empty()) {
        wr=initiator_wr_pool.front();
        initiator_wr_pool.pop_front();
    }
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return(wr);
}
static void wr_pool_target_push(gni_work_request *wr)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;

    log_debug(nnti_debug_level, "enter");

    memset(&wr->wc, 0, sizeof(nnti_gni_work_completion));
    wr->last_op       =0;
    wr->is_op_complete=FALSE;
    wr->op_state      =BUFFER_INIT;

    if (nthread_lock(&nnti_wr_pool_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    target_wr_pool.push_front(wr);
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return;
}
static void wr_pool_initiator_push(gni_work_request *wr)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;

    log_debug(nnti_debug_level, "enter");

    memset(&wr->wc, 0, sizeof(nnti_gni_work_completion));
    wr->last_op       =0;
    wr->is_op_complete=FALSE;
    wr->op_state      =BUFFER_INIT;

    if (nthread_lock(&nnti_wr_pool_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    initiator_wr_pool.push_front(wr);
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return;
}
static NNTI_result_t wr_pool_fini(void)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;
    gni_work_request *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    if (nthread_lock(&nnti_wr_pool_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    while (!target_wr_pool.empty()) {
        wr=target_wr_pool.front();
        target_wr_pool.pop_front();
        assert(wr);
        rc=wr_pool_deregister(wr);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to deregister target work request: rc=%d", rc);
            goto cleanup;
        }
        free(wr);
    }
    while (!initiator_wr_pool.empty()) {
        wr=initiator_wr_pool.front();
        initiator_wr_pool.pop_front();
        assert(wr);
        rc=wr_pool_deregister(wr);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to deregister initiator work request: rc=%d", rc);
            goto cleanup;
        }
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

    log_debug(debug_level, "enter (%d instance connections, %d peer connections)",
            connections_by_instance.size(), connections_by_peer.size());

    if (nthread_lock(&nnti_conn_instance_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    conn_by_inst_iter_t inst_iter;
    for (inst_iter = connections_by_instance.begin(); inst_iter != connections_by_instance.end(); inst_iter++) {
        log_debug(debug_level, "close connection (instance=%llu)", inst_iter->first);
        close_connection(inst_iter->second);
        connections_by_instance.erase(inst_iter);
    }
    nthread_unlock(&nnti_conn_instance_lock);

    if (nthread_lock(&nnti_conn_peer_lock)) log_warn(nnti_debug_level, "failed to get thread lock");
    conn_by_peer_iter_t peer_iter;
    for (peer_iter = connections_by_peer.begin(); peer_iter != connections_by_peer.end(); peer_iter++) {
        log_debug(debug_level, "close connection (peer.addr=%llu)", peer_iter->first.addr);
//        close_connection(peer_iter->second);
        connections_by_peer.erase(peer_iter);
    }
    nthread_unlock(&nnti_conn_peer_lock);

    log_debug(debug_level, "exit (%d instance connections, %d peer connections)",
            connections_by_instance.size(), connections_by_peer.size());

    return;
}

/**
 * @brief initialize
 */
static NNTI_result_t init_connection(
        gni_connection **conn,
        const int sock,
        const int is_server)
{
    int rc=NNTI_OK; /* return code */

    trios_declare_timer(call_time);

    log_debug(nnti_debug_level, "initializing gni connection");

    trios_start_timer(call_time);
    if (is_server) {
        rc = new_server_connection(*conn, sock);
    } else {
        rc = new_client_connection(*conn, sock);
    }
    trios_stop_timer("new connection", call_time);
    if (rc) {
        close_connection(*conn);
        goto out;
    }

    print_gni_conn(*conn);

out:
    return((NNTI_result_t)rc);
}

/*
 * At an explicit BYE message, or at finalize time, shut down a connection.
 * If descriptors are posted, defer and clean up the connection structures
 * later.
 */
static void close_connection(gni_connection *c)
{
    log_level debug_level = nnti_debug_level;  // nnti_ee_debug_level;

    if (c==NULL) return;

    log_debug(debug_level, "close_connection: start");

    print_gni_conn(c);

    if (c->peer_name) {
        free(c->peer_name);
        c->peer_name = NULL;
    }

    if (c->connection_type == CLIENT_CONNECTION) {
        client_req_queue_destroy(c);
    }

    c->state=DISCONNECTED;

    log_debug(debug_level, "close_connection: exit");
}

static int check_for_waiting_connection()
{
    NNTI_result_t rc = NNTI_OK;

    struct sockaddr_in ssin;
    socklen_t len;
    int s;
    gni_connection *conn = NULL;
    NNTI_peer_t peer;

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
        conn = (gni_connection *)calloc(1, sizeof(gni_connection));
        log_debug(nnti_debug_level, "calloc returned conn=%p.", conn);
        if (conn == NULL) {
            log_error(nnti_debug_level, "calloc returned NULL.  out of memory?: %s", strerror(errno));
            rc=NNTI_ENOMEM;
            goto cleanup;
        }

//        nthread_lock(&nnti_gni_lock);
        rc=init_connection(&conn, s, 1);
        if (rc!=NNTI_OK) {
            goto cleanup;
        }
        create_peer(
                &peer,
                conn->peer_name,
                conn->peer_addr,
                conn->peer_port,
                conn->peer_ptag,
                conn->peer_cookie,
                conn->peer_instance);

        insert_conn_peer(&peer, conn);
        insert_conn_instance(conn->peer_instance, conn);

        transition_connection_to_ready(s, conn);
//        nthread_unlock(&nnti_gni_lock);

        log_debug(nnti_debug_level, "accepted new connection from %s:%u", conn->peer_name, conn->peer_port);

        if (close(s) < 0) {
            log_error(nnti_debug_level, "failed to close new tcp socket");
        }

        if (logging_debug(nnti_debug_level)) {
            fprint_NNTI_peer(logger_get_file(), "peer",
                    "end of check_listen_socket_for_new_connections", &peer);
        }
    }

cleanup:
    if (rc != NNTI_OK) {
        if (conn!=NULL) free(conn);
    }
    return rc;
}

/**
 * Check for new connections.  The listening socket is left nonblocking
 * so this test can be quick; but accept is not really that quick compared
 * to polling an Gemini interface, for instance.  Returns >0 if an accept worked.
 */
static int check_listen_socket_for_new_connections()
{
    bool done=false;
    while(!done) {
        if (check_for_waiting_connection() != NNTI_OK) {
            done=true;
        }
    }

    return(NNTI_OK);
}

static uint32_t get_cpunum(void)
{
  int i, j;
  uint32_t cpu_num;

  cpu_set_t coremask;

  (void)sched_getaffinity(0, sizeof(coremask), &coremask);

  for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, &coremask)) {
      int run = 0;
      for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, &coremask)) run++;
        else break;
      }
      if (!run) {
        cpu_num=i;
      } else {
        fprintf(stdout, "This thread is bound to multiple CPUs(%d).  Using lowest numbered CPU(%d).", run+1, cpu_num);
        cpu_num=i;
      }
    }
  }
  return(cpu_num);
}


static void get_alps_info(
        alpsAppGni_t *alps_info)
{
    int alps_rc=0;
    int req_rc=0;
    size_t rep_size=0;

    uint64_t apid=0;
    alpsAppLLIGni_t *alps_info_list;
    char buf[1024];

    alps_info_list=(alpsAppLLIGni_t *)&buf[0];

    alps_app_lli_lock();

    log_debug(nnti_debug_level, "sending ALPS request");
    alps_rc = alps_app_lli_put_request(ALPS_APP_LLI_ALPS_REQ_GNI, NULL, 0);
    if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_put_request failed: %d", alps_rc);
    log_debug(nnti_debug_level, "waiting for ALPS reply");
    alps_rc = alps_app_lli_get_response(&req_rc, &rep_size);
    if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response failed: alps_rc=%d", alps_rc);
    if (req_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response failed: req_rc=%d", req_rc);
    if (rep_size != 0) {
        log_debug(nnti_debug_level, "waiting for ALPS reply bytes (%d) ; sizeof(alps_info)==%d ; sizeof(alps_info_list)==%d", rep_size, sizeof(alps_info), sizeof(alps_info_list));
        alps_rc = alps_app_lli_get_response_bytes(alps_info_list, rep_size);
        if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response_bytes failed: %d", alps_rc);
    }

    log_debug(nnti_debug_level, "sending ALPS request");
    alps_rc = alps_app_lli_put_request(ALPS_APP_LLI_ALPS_REQ_APID, NULL, 0);
    if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_put_request failed: %d", alps_rc);
    log_debug(nnti_debug_level, "waiting for ALPS reply");
    alps_rc = alps_app_lli_get_response(&req_rc, &rep_size);
    if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response failed: alps_rc=%d", alps_rc);
    if (req_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response failed: req_rc=%d", req_rc);
    if (rep_size != 0) {
        log_debug(nnti_debug_level, "waiting for ALPS reply bytes (%d) ; sizeof(apid)==%d", rep_size, sizeof(apid));
        alps_rc = alps_app_lli_get_response_bytes(&apid, rep_size);
        if (alps_rc != 0) log_debug(nnti_debug_level, "alps_app_lli_get_response_bytes failed: %d", alps_rc);
    }

    alps_app_lli_unlock();

    memcpy(alps_info, (alpsAppGni_t*)alps_info_list->u.buf, sizeof(alpsAppGni_t));
    transport_global_data.apid=apid;

    log_debug(nnti_debug_level, "apid                 =%llu", (unsigned long long)apid);
    log_debug(nnti_debug_level, "alps_info->device_id =%llu", (unsigned long long)alps_info->device_id);
    log_debug(nnti_debug_level, "alps_info->local_addr=%lld", (long long)alps_info->local_addr);
    log_debug(nnti_debug_level, "alps_info->cookie    =%llu", (unsigned long long)alps_info->cookie);
    log_debug(nnti_debug_level, "alps_info->ptag      =%llu", (unsigned long long)alps_info->ptag);

    log_debug(nnti_debug_level, "ALPS response - apid(%llu) alps_info->device_id(%llu) alps_info->local_addr(%llu) "
            "alps_info->cookie(%llu) alps_info->ptag(%llu)",
            (unsigned long long)apid,
            (unsigned long long)alps_info->device_id,
            (long long)alps_info->local_addr,
            (unsigned long long)alps_info->cookie,
            (unsigned long long)alps_info->ptag);

    return;
}

static void print_wc(const nnti_gni_work_completion *wc)
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    log_debug(nnti_debug_level, "wc=%p, wc.op=%d, wc.inst_id=%llu, wc.byte_len=%llu, wc.byte_offset=%llu, wc.src_offset=%llu, wc.dest_offset=%llu",
            wc,
            wc->op,
            wc->inst_id,
            wc->byte_len,
            wc->byte_offset,
            wc->src_offset,
            wc->dest_offset);
}

static void print_cq_event(
        const gni_cq_entry_t *event)
{
    if (gni_cq_get_status(*event) != 0) {
        log_error(nnti_debug_level, "event=%p, event.data=%llu, event.source=%llu, event.status=%llu, "
                "event.info=%llu, event.overrun=%llu, event.inst_id=%llu, event.tid=%llu, event.msg_id=%llu, event.type=%llu",
                event,
                (uint64_t)gni_cq_get_data(*event),
                (uint64_t)gni_cq_get_source(*event),
                (uint64_t)gni_cq_get_status(*event),
                (uint64_t)gni_cq_get_info(*event),
                (uint64_t)gni_cq_overrun(*event),
                (uint64_t)gni_cq_get_inst_id(*event),
                (uint64_t)gni_cq_get_tid(*event),
                (uint64_t)gni_cq_get_msg_id(*event),
                (uint64_t)gni_cq_get_type(*event));
    } else {
        log_debug(nnti_debug_level, "event=%p, event.data=%llu, event.source=%llu, event.status=%llu, "
                "event.info=%llu, event.overrun=%llu, event.inst_id=%llu, event.tid=%llu, event.msg_id=%llu, event.type=%llu",
                event,
                (uint64_t)gni_cq_get_data(*event),
                (uint64_t)gni_cq_get_source(*event),
                (uint64_t)gni_cq_get_status(*event),
                (uint64_t)gni_cq_get_info(*event),
                (uint64_t)gni_cq_overrun(*event),
                (uint64_t)gni_cq_get_inst_id(*event),
                (uint64_t)gni_cq_get_tid(*event),
                (uint64_t)gni_cq_get_msg_id(*event),
                (uint64_t)gni_cq_get_type(*event));
    }
}

static void print_post_desc(
        const gni_post_descriptor_t *post_desc_ptr)
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    if (post_desc_ptr != NULL) {
        log_debug(nnti_debug_level, "post_desc_ptr                  ==%p", (uint64_t)post_desc_ptr);
        log_debug(nnti_debug_level, "post_desc_ptr->next_descr      ==%p", (uint64_t)post_desc_ptr->next_descr);
        log_debug(nnti_debug_level, "post_desc_ptr->prev_descr      ==%p", (uint64_t)post_desc_ptr->prev_descr);

        log_debug(nnti_debug_level, "post_desc_ptr->post_id         ==%llu", (uint64_t)post_desc_ptr->post_id);
        log_debug(nnti_debug_level, "post_desc_ptr->status          ==%llu", (uint64_t)post_desc_ptr->status);
        log_debug(nnti_debug_level, "post_desc_ptr->cq_mode_complete==%llu", (uint64_t)post_desc_ptr->cq_mode_complete);

        log_debug(nnti_debug_level, "post_desc_ptr->type            ==%llu", (uint64_t)post_desc_ptr->type);
        log_debug(nnti_debug_level, "post_desc_ptr->cq_mode         ==%llu", (uint64_t)post_desc_ptr->cq_mode);
        log_debug(nnti_debug_level, "post_desc_ptr->dlvr_mode       ==%llu", (uint64_t)post_desc_ptr->dlvr_mode);
        log_debug(nnti_debug_level, "post_desc_ptr->local_addr      ==%llu", (uint64_t)post_desc_ptr->local_addr);
        log_debug(nnti_debug_level, "post_desc_ptr->remote_addr     ==%llu", (uint64_t)post_desc_ptr->remote_addr);
        log_debug(nnti_debug_level, "post_desc_ptr->length          ==%llu", (uint64_t)post_desc_ptr->length);
        log_debug(nnti_debug_level, "post_desc_ptr->rdma_mode       ==%llu", (uint64_t)post_desc_ptr->rdma_mode);
        log_debug(nnti_debug_level, "post_desc_ptr->src_cq_hndl     ==%llu", (uint64_t)post_desc_ptr->src_cq_hndl);
        log_debug(nnti_debug_level, "post_desc_ptr->sync_flag_value ==%llu", (uint64_t)post_desc_ptr->sync_flag_value);
        log_debug(nnti_debug_level, "post_desc_ptr->sync_flag_addr  ==%llu", (uint64_t)post_desc_ptr->sync_flag_addr);
        log_debug(nnti_debug_level, "post_desc_ptr->amo_cmd         ==%llu", (uint64_t)post_desc_ptr->amo_cmd);
        log_debug(nnti_debug_level, "post_desc_ptr->first_operand   ==%llu", (uint64_t)post_desc_ptr->first_operand);
        log_debug(nnti_debug_level, "post_desc_ptr->second_operand  ==%llu", (uint64_t)post_desc_ptr->second_operand);
        log_debug(nnti_debug_level, "post_desc_ptr->cqwrite_value   ==%llu", (uint64_t)post_desc_ptr->cqwrite_value);
    } else {
        log_debug(nnti_debug_level, "post_desc_ptr == NULL");
    }
}

static void print_gni_conn(gni_connection *c)
{
    log_level debug_level=nnti_debug_level;

    if (!logging_debug(debug_level)) {
        return;
    }

    log_debug(debug_level, "c->peer_name       =%s", c->peer_name);
    log_debug(debug_level, "c->peer_addr       =%u", c->peer_addr);
    log_debug(debug_level, "c->peer_port       =%u", (uint32_t)c->peer_port);
    log_debug(debug_level, "c->peer_instance   =%llu", (uint64_t)c->peer_instance);

    log_debug(debug_level, "c->state           =%d", c->state);
}

//static void print_put_buf(void *buf, uint32_t size)
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
//        log_debug(nnti_select_debug_level, "array[%d].int_val=%u ; array[%d].float_val=%f ; array[%d].double_val=%f",
//                idx, array[idx].int_val,
//                idx, array[idx].float_val,
//                idx, array[idx].double_val);
//    }
//
//}

static void print_raw_buf(void *buf, uint32_t size)
{
    if (logging_debug(nnti_debug_level)) {
        FILE* f=logger_get_file();
        u_int64_t print_limit=(size<90) ? size : 90;
        fprintf(f, "\nbuf (%p)\n", buf);
        fflush(f);
        if (buf != NULL) {
            int l=0;
            for (l=0;l<print_limit;l++) {
                if (l%30 == 0) fprintf(f, "\nbuf (%lu) (offset(%d)) => ", buf, l);
                fprintf(f, "%02hhX", ((char *)buf)[l]);
            }
            fprintf(f, "\n");
        }
    }
}

static uint16_t get_dlvr_mode_from_env()
{
    char *mode=getenv("NSSI_GNI_DELIVERY_MODE");

    log_debug(nnti_debug_level, "NSSI_GNI_DELIVERY_MODE=%s", mode);
    if ((mode==NULL) || !strcmp(mode, "GNI_DLVMODE_PERFORMANCE")) {
        log_debug(nnti_debug_level, "setting delivery mode to GNI_DLVMODE_PERFORMANCE");
        return GNI_DLVMODE_PERFORMANCE;
    } else if (!strcmp(mode, "GNI_DLVMODE_IN_ORDER")) {
        log_debug(nnti_debug_level, "setting delivery mode to GNI_DLVMODE_IN_ORDER");
        return GNI_DLVMODE_IN_ORDER;
    } else if (!strcmp(mode, "GNI_DLVMODE_NO_ADAPT")) {
        log_debug(nnti_debug_level, "setting delivery mode to GNI_DLVMODE_NO_ADAPT");
        return GNI_DLVMODE_NO_ADAPT;
    } else if (!strcmp(mode, "GNI_DLVMODE_NO_HASH")) {
        log_debug(nnti_debug_level, "setting delivery mode to GNI_DLVMODE_NO_HASH");
        return GNI_DLVMODE_NO_HASH;
    } else if (!strcmp(mode, "GNI_DLVMODE_NO_RADAPT")) {
        log_debug(nnti_debug_level, "setting delivery mode to GNI_DLVMODE_NO_RADAPT");
        return GNI_DLVMODE_NO_RADAPT;
    } else {
        log_debug(nnti_debug_level, "defaulting delivery mode to GNI_DLVMODE_PERFORMANCE");
        return GNI_DLVMODE_PERFORMANCE;
    }
}
static void set_dlvr_mode(
        gni_post_descriptor_t *pd)
{
    pd->dlvr_mode=transport_global_data.delivery_mode;
}

static void set_rdma_mode(
        gni_post_descriptor_t *pd)
{
    if (config.use_rdma_fence==true) {
        pd->rdma_mode=GNI_RDMAMODE_FENCE;
    }
}

static void set_post_desc(gni_post_descriptor_t *pd, gni_buffer_type target_buf_type, uint32_t buf_length)
{
    if (target_buf_type == GET_SRC_BUFFER) {
        // set type
        switch (config.rdma_mode) {
            case RDMA_FMA:
                pd->type=GNI_POST_FMA_GET;
                break;
            case RDMA_CROSSOVER:
                if (buf_length > config.fma_bte_crossover_size) {
                    pd->type=GNI_POST_RDMA_GET;
                } else {
                    pd->type=GNI_POST_FMA_GET;
                }
                break;
            case RDMA_BTE:
            case RDMA_MIXED:
            default:
                pd->type=GNI_POST_RDMA_GET;
                break;
        }
        // set cq_mode
        if (config.use_rdma_events==true) {
            pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT|GNI_CQMODE_REMOTE_EVENT;
        } else {
            pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT;
        }

    } else if (target_buf_type == PUT_DST_BUFFER) {
        // set type
        switch (config.rdma_mode) {
            case RDMA_FMA:
                pd->type=GNI_POST_FMA_PUT;
                break;
            case RDMA_CROSSOVER:
                if (buf_length > config.fma_bte_crossover_size) {
                    pd->type=GNI_POST_RDMA_PUT;
                } else {
                    pd->type=GNI_POST_FMA_PUT;
                }
                break;
            case RDMA_BTE:
            case RDMA_MIXED:
            default:
                pd->type=GNI_POST_RDMA_PUT;
                break;
        }
        // set cq_mode
        if (config.use_rdma_events==true) {
            pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT|GNI_CQMODE_REMOTE_EVENT;
        } else {
            pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT;
        }

    } else if (target_buf_type == RECEIVE_BUFFER) {
        // set type
        switch (config.rdma_mode) {
            case RDMA_CROSSOVER:
                if (buf_length > config.fma_bte_crossover_size) {
                    pd->type=GNI_POST_RDMA_PUT;
                } else {
                    pd->type=GNI_POST_FMA_PUT;
                }
                break;
            case RDMA_BTE:
                pd->type=GNI_POST_RDMA_PUT;
                break;
            case RDMA_FMA:
            case RDMA_MIXED:
            default:
                pd->type=GNI_POST_FMA_PUT;
                break;
        }
        // set cq_mode
        pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT|GNI_CQMODE_REMOTE_EVENT;

    } else if (target_buf_type == REQUEST_BUFFER) {
        // set type
        switch (config.rdma_mode) {
            case RDMA_CROSSOVER:
                if (buf_length > config.fma_bte_crossover_size) {
                    pd->type=GNI_POST_RDMA_PUT;
                } else {
                    pd->type=GNI_POST_FMA_PUT;
                }
                break;
            case RDMA_BTE:
                pd->type=GNI_POST_RDMA_PUT;
                break;
            case RDMA_FMA:
            case RDMA_MIXED:
            default:
                pd->type=GNI_POST_FMA_PUT;
                break;
        }
        // set cq_mode
        pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT;
    } else {
        log_error(nnti_debug_level, "unknown target buffer type(%llu).", (uint64_t)target_buf_type);
    }

    set_dlvr_mode(pd);
    set_rdma_mode(pd);
}

static void set_wc_post_desc(gni_post_descriptor_t *pd, uint32_t buf_length)
{
    // set type
    switch (config.rdma_mode) {
        case RDMA_CROSSOVER:
            if (buf_length > config.fma_bte_crossover_size) {
                pd->type=GNI_POST_RDMA_PUT;
            } else {
                pd->type=GNI_POST_FMA_PUT;
            }
            break;
        case RDMA_BTE:
            pd->type=GNI_POST_RDMA_PUT;
            break;
        case RDMA_FMA:
        case RDMA_MIXED:
        default:
            pd->type=GNI_POST_FMA_PUT;
            break;
    }
    // set cq_mode
    pd->cq_mode=GNI_CQMODE_GLOBAL_EVENT|GNI_CQMODE_REMOTE_EVENT;

    set_dlvr_mode(pd);
    set_rdma_mode(pd);
}

static int post_wait(
        gni_cq_handle_t cq_hdl,
        int             timeout,
        int             retries)
{
    int rc=0;
    int i=0;
    trios_declare_timer(call_time);
    gni_post_descriptor_t *post_desc_ptr;
    gni_cq_entry_t ev_data;

    log_debug(nnti_ee_debug_level, "enter");

    memset(&ev_data, 0, sizeof(ev_data));
    for(i=0;i<=retries;i++) {
        log_debug(nnti_debug_level, "calling CqWaitEvent");
        trios_start_timer(call_time);
        rc=GNI_CqWaitEvent_wrapper(cq_hdl, timeout, &ev_data);
        trios_stop_timer("post_wait - CqWaitEvent", call_time);
        if (rc==GNI_RC_SUCCESS) {
            break;
        } else if (rc!=GNI_RC_TIMEOUT) {
            log_error(nnti_debug_level, "CqWaitEvent failed: %d", rc);
        }
    }

    log_debug(nnti_debug_level, "calling GetComplete");
    trios_start_timer(call_time);
    rc=GNI_GetCompleted (cq_hdl, ev_data, &post_desc_ptr);
    trios_stop_timer("post_wait - GetCompleted", call_time);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "GetCompleted failed: %d", rc);
    print_post_desc(post_desc_ptr);

    log_debug(nnti_ee_debug_level, "exit");

    return(rc);
}


static int reset_req_index(
        gni_request_queue_handle  *req_queue_attrs)
{
    int rc=0;
    gni_post_descriptor_t  post_desc;

    uint64_t value_before_reset=0;
    uint64_t value_before_reset_addr=(uint64_t)&value_before_reset;
    gni_mem_handle_t value_before_reset_mem_hdl;
    gni_cq_handle_t reset_cq_hdl;
    gni_ep_handle_t reset_ep_hdl;

    log_debug(nnti_ee_debug_level, "enter");

    rc=GNI_CqCreate (transport_global_data.nic_hdl, 1, 0, GNI_CQ_BLOCKING, NULL, NULL, &reset_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqCreate(value_before_reset_cq_hdl) failed: %d", rc);

    rc=GNI_MemRegister (transport_global_data.nic_hdl, value_before_reset_addr, sizeof(uint64_t), NULL, GNI_MEM_READWRITE, (uint32_t)-1, &value_before_reset_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemRegister(value_before_reset) failed: %d", rc);

    rc=GNI_EpCreate (transport_global_data.nic_hdl, reset_cq_hdl, &reset_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpCreate(reset_ep_hdl) failed: %d", rc);
    rc=GNI_EpBind (reset_ep_hdl, transport_global_data.alps_info.local_addr, transport_global_data.instance);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpBind(reset_ep_hdl) failed: %d", rc);

//    log_debug(nnti_debug_level, "index before reset(%llu).", req_queue_attrs->req_index);


    memset(&post_desc, 0, sizeof(gni_post_descriptor_t));
    post_desc.type           =GNI_POST_AMO;
    post_desc.cq_mode        =GNI_CQMODE_GLOBAL_EVENT;

    set_dlvr_mode(&post_desc);
    set_rdma_mode(&post_desc);

    post_desc.local_addr     =value_before_reset_addr;
    post_desc.local_mem_hndl =value_before_reset_mem_hdl;
    post_desc.remote_addr    =req_queue_attrs->req_index_addr;
    post_desc.remote_mem_hndl=req_queue_attrs->req_index_mem_hdl;
    post_desc.length         =sizeof(uint64_t);
    post_desc.amo_cmd        =GNI_FMA_ATOMIC_FAND;
    post_desc.first_operand  =0;

    log_debug(nnti_debug_level, "calling PostFma(reset index ep_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
            reset_ep_hdl, post_desc.local_addr, post_desc.remote_addr);
    rc=GNI_PostFma(reset_ep_hdl, &post_desc);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostFma(reset index) failed: %d", rc);

    post_wait(reset_cq_hdl, 1000, 0);

    rc=GNI_MemDeregister (transport_global_data.nic_hdl, &value_before_reset_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemDeregister(1) failed: %d", rc);

    rc=GNI_EpUnbind (reset_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpUnbind() failed: %d", rc);
    rc=GNI_EpDestroy (reset_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpDestroy() failed: %d", rc);

    rc=GNI_CqDestroy (reset_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqDestroy() failed: %d", rc);

    req_queue_attrs->last_index_before_reset=value_before_reset;
    req_queue_attrs->req_processed_reset_limit=value_before_reset;
    log_debug(nnti_event_debug_level, "index before reset(%llu).", req_queue_attrs->last_index_before_reset);
    log_debug(nnti_event_debug_level, "index after reset(%llu).", req_queue_attrs->req_index);

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int fetch_add_buffer_offset(
        nnti_gni_client_queue       *local_req_queue_attrs,
        nnti_gni_recv_queue_attrs *remote_req_queue_attrs,
        uint64_t                     addend,
        uint64_t                    *prev_offset)
{
    int rc=0;
    trios_declare_timer(call_time);
    gni_post_descriptor_t  post_desc;
    gni_cq_entry_t ev_data;
    int extras_absorbed=0;

    log_level debug_level=nnti_debug_level;

//    static uint64_t last_offset=0;

    log_debug(nnti_ee_debug_level, "enter");

    memset(&post_desc, 0, sizeof(gni_post_descriptor_t));
    post_desc.type           =GNI_POST_AMO;
    post_desc.cq_mode        =GNI_CQMODE_GLOBAL_EVENT;

    set_dlvr_mode(&post_desc);
    set_rdma_mode(&post_desc);

    post_desc.local_addr     =local_req_queue_attrs->req_index_addr;
    post_desc.local_mem_hndl =local_req_queue_attrs->req_index_mem_hdl;
    post_desc.remote_addr    =remote_req_queue_attrs->req_index_addr;
    post_desc.remote_mem_hndl=remote_req_queue_attrs->req_index_mem_hdl;
    post_desc.length         =sizeof(uint64_t);
    post_desc.amo_cmd        =GNI_FMA_ATOMIC_FADD;
    post_desc.first_operand  =addend;
    post_desc.second_operand =0;

    do {
        log_debug(debug_level, "calling PostFma(fetch add - req_index_ep_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                local_req_queue_attrs->req_index_ep_hdl, post_desc.local_addr, (uint64_t)post_desc.remote_addr);
        rc=GNI_PostFma(local_req_queue_attrs->req_index_ep_hdl, &post_desc);
        if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "PostFma(fetch add) failed: %d", rc);

        post_wait(local_req_queue_attrs->req_index_cq_hdl, 1000, 0);

        log_debug(debug_level, "fetched queue_index(%llu)", local_req_queue_attrs->req_index);
        *prev_offset=local_req_queue_attrs->req_index;
        if (*prev_offset >= remote_req_queue_attrs->req_count) {
            uint64_t  ui64;
            uint32_t *ptr32=NULL;


            log_debug(debug_level, "fetched queue_index(%llu) >= req_count(%llu)",
                    local_req_queue_attrs->req_index, remote_req_queue_attrs->req_count);

            memset(&ev_data, 0, sizeof(ev_data));
            do {
                log_debug(debug_level, "calling CqWaitEvent(unblock)");
                trios_start_timer(call_time);
                rc=GNI_CqWaitEvent_wrapper (local_req_queue_attrs->unblock_mem_cq_hdl, 1000, &ev_data);
                trios_stop_timer("unblock", call_time);
                if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "CqWaitEvent(unblock) failed: %d", rc);
            } while (rc!=GNI_RC_SUCCESS);
            ui64=gni_cq_get_data(ev_data);
            ptr32=(uint32_t*)&ui64;
            ptr32++;
            log_debug(debug_level, "received unblock from server (ACKing requests <= %llu)", *ptr32);

            if (rc==GNI_RC_SUCCESS) {
                extras_absorbed=0;
                log_debug(debug_level, "calling CqGetEvent(absorb extra unblocks)");
                do {
                    log_debug(debug_level, "calling CqGetEvent(absorb extra unblocks)");
                    rc=GNI_CqGetEvent (local_req_queue_attrs->unblock_mem_cq_hdl, &ev_data);
                    if (rc==GNI_RC_SUCCESS) {
                        extras_absorbed++;

                        ui64=gni_cq_get_data(ev_data);
                        ptr32=(uint32_t*)&ui64;
                        ptr32++;
                        log_debug(debug_level, "received extra unblock from server (ACKing requests <= %llu)", *ptr32);
                    }
                } while (rc==GNI_RC_SUCCESS);
                if ((rc!=GNI_RC_SUCCESS) && (rc!=GNI_RC_NOT_DONE)) log_error(debug_level, "CqGetEvent(absorb extra unblocks) failed: %d", rc);
                log_debug(debug_level, "absorbed %d extra unblocks)", extras_absorbed);
            }
        } else if (*prev_offset < local_req_queue_attrs->last_offset) {
            uint64_t  ui64;
            uint32_t *ptr32=NULL;

            /* absorb one unblock message */
            log_debug(nnti_event_debug_level, "calling CqGetEvent(absorb one unblock)");
            rc=GNI_CqWaitEvent_wrapper (local_req_queue_attrs->unblock_mem_cq_hdl, 1000, &ev_data);
            if (rc==GNI_RC_SUCCESS) {
                ui64=gni_cq_get_data(ev_data);
                ptr32=(uint32_t*)&ui64;
                ptr32++;
                log_debug(nnti_event_debug_level, "received one unblock from server (ACKing requests <= %llu)", *ptr32);
            } else if (rc!=GNI_RC_NOT_DONE) {
                log_error(debug_level, "CqGetEvent(absorb one unblock) failed: %d", rc);
            }
        }
    } while (*prev_offset >= remote_req_queue_attrs->req_count);

    local_req_queue_attrs->last_offset=*prev_offset;

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int send_req(
        const NNTI_peer_t           *peer_hdl,
        nnti_gni_client_queue       *local_req_queue_attrs,
        nnti_gni_recv_queue_attrs *remote_req_queue_attrs,
        uint64_t                     offset,
        const NNTI_buffer_t         *reg_buf,
        gni_work_request            *wr)
{
    int rc=0;
    gni_memory_handle *gni_mem_hdl=NULL;

    trios_declare_timer(call_time);

    bool use_fma=true;

    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    memset(&wr->post_desc, 0, sizeof(gni_post_descriptor_t));
    set_post_desc(&wr->post_desc, REQUEST_BUFFER, reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.size);

    wr->last_op =GNI_OP_SEND_REQUEST;
    wr->op_state=RDMA_WRITE_INIT;

    wr->post_desc.local_addr           =reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf;
    wr->post_desc.local_mem_hndl.qword1=reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
    wr->post_desc.local_mem_hndl.qword2=reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword2;
    wr->post_desc.remote_addr          =remote_req_queue_attrs->req_buffer_addr+offset;
    wr->post_desc.remote_mem_hndl      =remote_req_queue_attrs->req_mem_hdl;
    wr->post_desc.length               =reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.size;

    print_raw_buf((void *)wr->post_desc.local_addr, wr->post_desc.length);

    GNI_EpSetEventData(
            local_req_queue_attrs->req_ep_hdl,
            hash6432shift((uint64_t)reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf),
            offset);

    if (config.rdma_mode==RDMA_CROSSOVER) {
        if (wr->post_desc.length > config.fma_bte_crossover_size) {
            use_fma=false;
        } else {
            use_fma=true;
        }
    } else if (config.rdma_mode==RDMA_BTE) {
        use_fma=false;
    } else if ((config.rdma_mode==RDMA_FMA) || (config.rdma_mode==RDMA_MIXED)) {
        use_fma=true;
    }

    if (use_fma) {
        log_debug(nnti_debug_level, "calling PostFma(send req ep_hdl(%llu), cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                local_req_queue_attrs->req_ep_hdl, local_req_queue_attrs->req_cq_hdl, wr->post_desc.local_addr, wr->post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostFma(local_req_queue_attrs->req_ep_hdl, &wr->post_desc);
        trios_stop_timer("PostFma req", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostFma(fma put) failed: %d", rc);
    } else {
        log_debug(nnti_debug_level, "calling PostRdma(send req ep_hdl(%llu), cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                local_req_queue_attrs->req_ep_hdl, local_req_queue_attrs->req_cq_hdl, wr->post_desc.local_addr, wr->post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostRdma(local_req_queue_attrs->req_ep_hdl, &wr->post_desc);
        trios_stop_timer("PostRdma req", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostRdma(rdma put) failed: %d", rc);
    }

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int send_req_wc(
        const NNTI_peer_t           *peer_hdl,
        nnti_gni_client_queue       *local_req_queue_attrs,
        nnti_gni_recv_queue_attrs *remote_req_queue_attrs,
        uint64_t                     offset,
        const NNTI_buffer_t         *reg_buf,
        gni_work_request            *wr)
{
    int rc=0;
    gni_ep_handle_t    ep_hdl;
    gni_memory_handle *gni_mem_hdl=NULL;

    bool use_fma=true;

    trios_declare_timer(call_time);

    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
    assert(gni_mem_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    print_wc(&wr->wc);

    memset(&wr->wc_post_desc, 0, sizeof(gni_post_descriptor_t));
    set_wc_post_desc(&wr->wc_post_desc, sizeof(nnti_gni_work_completion));

    wr->wc_post_desc.local_addr     =(uint64_t)&wr->wc;
    wr->wc_post_desc.local_mem_hndl =wr->wc_mem_hdl;
    wr->wc_post_desc.remote_addr    =remote_req_queue_attrs->wc_buffer_addr+offset;
    wr->wc_post_desc.remote_mem_hndl=remote_req_queue_attrs->wc_mem_hdl;
    wr->wc_post_desc.length         =sizeof(nnti_gni_work_completion);

    print_raw_buf((void *)wr->wc_post_desc.local_addr, wr->wc_post_desc.length);

    gni_connection *conn=get_conn_peer(peer_hdl);
    assert(conn);
    ep_hdl=conn->ep_hdl;

    wr->peer_instance=conn->peer_instance;

    GNI_EpSetEventData(
            ep_hdl,
            hash6432shift((uint64_t)reg_buf->buffer_addr.NNTI_remote_addr_t_u.gni.buf),
            offset/sizeof(nnti_gni_work_completion));

    if (config.rdma_mode==RDMA_CROSSOVER) {
        if (wr->wc_post_desc.length > config.fma_bte_crossover_size) {
            use_fma=false;
        } else {
            use_fma=true;
        }
    } else if (config.rdma_mode==RDMA_BTE) {
        use_fma=false;
    } else if ((config.rdma_mode==RDMA_FMA) || (config.rdma_mode==RDMA_MIXED)) {
        use_fma=true;
    }

    if (use_fma) {
        log_debug(nnti_debug_level, "calling PostFma(send_req_wc ep_hdl(%llu) wc_cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                ep_hdl, transport_global_data.ep_cq_hdl, wr->wc_post_desc.local_addr, wr->wc_post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostFma(ep_hdl, &wr->wc_post_desc);
        trios_stop_timer("PostFma wc", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostFma(fma put) failed: %d", rc);
    } else {
        log_debug(nnti_debug_level, "calling PostRdma(send_req_wc ep_hdl(%llu) wc_cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                ep_hdl, transport_global_data.ep_cq_hdl, wr->wc_post_desc.local_addr, wr->wc_post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostRdma(ep_hdl, &wr->wc_post_desc);
        trios_stop_timer("PostRdma wc", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostRdma(rdma put) failed: %d", rc);
    }

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int request_send(
        const NNTI_peer_t           *peer_hdl,
        nnti_gni_client_queue       *client_q,
        nnti_gni_recv_queue_attrs *server_q,
        const NNTI_buffer_t         *reg_buf,
        int                          req_num)
{
    int rc=0;
    uint64_t offset=0;
    trios_declare_timer(call_time);

    gni_memory_handle *gni_mem_hdl=NULL;
    gni_work_request  *wr=NULL;

    uint32_t wc_size=sizeof(nnti_gni_work_completion);

    log_debug(nnti_ee_debug_level, "enter");

    gni_mem_hdl=(gni_memory_handle *)reg_buf->transport_private;
    assert(gni_mem_hdl);
    if (config.use_wr_pool) {
        wr=wr_pool_initiator_pop();
    } else {
        wr=(gni_work_request *)calloc(1, sizeof(gni_work_request));
    }
    assert(wr);

    wr->reg_buf = reg_buf;

    log_debug(nnti_debug_level, "calling fetch_add_buffer_offset()");
    trios_start_timer(call_time);
    rc=fetch_add_buffer_offset(client_q, server_q, 1, &offset);
    trios_stop_timer("fetch_add_buffer_offset", call_time);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "fetch_add_buffer_offset() failed: %d", rc);

    if (!config.use_wr_pool) {
        wr->wc_registered  =FALSE;
    }
    wr->wc.ack_received=0;
    wr->wc.inst_id     =transport_global_data.instance;
    wr->wc.op         =GNI_OP_SEND_REQUEST;
    wr->wc.byte_len   =reg_buf->payload_size;
    wr->wc.byte_offset=offset*server_q->req_size;
    wr->wc.src_offset =0;
    wr->wc.dest_offset=offset*server_q->req_size;

    nthread_lock(&gni_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "calling send_req()");
    trios_start_timer(call_time);
    rc=send_req(peer_hdl, client_q, server_q, offset*server_q->req_size, reg_buf, wr);
    trios_stop_timer("send_req", call_time);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "send_req() failed: %d", rc);

    log_debug(nnti_debug_level, "calling send_req_wc()");
    trios_start_timer(call_time);
    rc=send_req_wc(peer_hdl, client_q, server_q, offset*wc_size, reg_buf, wr);
    trios_stop_timer("send_req_wc", call_time);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "send_req_wc() failed: %d", rc);

    gni_mem_hdl->wr_queue.push_back(wr);
    insert_wr_wrhash(wr);

    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

    wr->last_op      =GNI_OP_SEND_REQUEST;
    wr->peer_instance=peer_hdl->peer.NNTI_remote_process_t_u.gni.inst_id;
    log_debug(nnti_event_debug_level, "wr->last_peer==%lu", wr->peer_instance);

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int send_buffer(
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *src_hdl,
        const NNTI_buffer_t *dest_hdl,
        gni_work_request    *wr)
{
    int rc=0;

    gni_ep_handle_t        ep_hdl;

    gni_memory_handle *gni_mem_hdl=NULL;

    bool use_fma=true;

    trios_declare_timer(call_time);

    gni_mem_hdl=(gni_memory_handle *)src_hdl->transport_private;
    assert(gni_mem_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    memset(&wr->post_desc, 0, sizeof(gni_post_descriptor_t));
    set_post_desc(&wr->post_desc, RECEIVE_BUFFER, src_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.size);

    wr->last_op=GNI_OP_SEND_BUFFER;
    wr->op_state=RDMA_WRITE_INIT;

    wr->post_desc.local_addr            =src_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf;
    wr->post_desc.local_mem_hndl.qword1 =src_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
    wr->post_desc.local_mem_hndl.qword2 =src_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword2;
    wr->post_desc.remote_addr           =dest_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf;
    wr->post_desc.remote_mem_hndl.qword1=dest_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword1;
    wr->post_desc.remote_mem_hndl.qword2=dest_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword2;
    wr->post_desc.length                =src_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.size;

    print_raw_buf((void *)wr->post_desc.local_addr, wr->post_desc.length);

    gni_connection *conn=get_conn_peer(peer_hdl);
    assert(conn);
    ep_hdl=conn->ep_hdl;

    GNI_EpSetEventData(
            ep_hdl,
            hash6432shift((uint64_t)src_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf),
            hash6432shift((uint64_t)dest_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf));

    if (config.rdma_mode==RDMA_CROSSOVER) {
        if (wr->post_desc.length > config.fma_bte_crossover_size) {
            use_fma=false;
        } else {
            use_fma=true;
        }
    } else if (config.rdma_mode==RDMA_BTE) {
        use_fma=false;
    } else if ((config.rdma_mode==RDMA_FMA) || (config.rdma_mode==RDMA_MIXED)) {
        use_fma=true;
    }

    if (use_fma) {
        log_debug(nnti_debug_level, "calling PostFma(send buffer ep_hdl(%llu), cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                ep_hdl, transport_global_data.ep_cq_hdl, wr->post_desc.local_addr, wr->post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostFma(ep_hdl, &wr->post_desc);
        trios_stop_timer("PostFma req", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostFma(fma put) failed: %d", rc);
    } else {
        log_debug(nnti_debug_level, "calling PostRdma(send buffer ep_hdl(%llu), cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                ep_hdl, transport_global_data.ep_cq_hdl, wr->post_desc.local_addr, wr->post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostRdma(ep_hdl, &wr->post_desc);
        trios_stop_timer("PostRdma req", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostRdma(rdma put) failed: %d", rc);
    }

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int send_buffer_wc(
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *src_hdl,
        const NNTI_buffer_t *dest_hdl,
        gni_work_request    *wr)
{
    int rc=0;
    gni_ep_handle_t    ep_hdl;
    gni_memory_handle *gni_mem_hdl=NULL;

    bool use_fma=true;

    trios_declare_timer(call_time);

    gni_mem_hdl=(gni_memory_handle *)src_hdl->transport_private;
    assert(gni_mem_hdl);

    log_debug(nnti_ee_debug_level, "enter");

    print_wc(&wr->wc);

    memset(&wr->wc_post_desc, 0, sizeof(gni_post_descriptor_t));
    set_wc_post_desc(&wr->wc_post_desc, sizeof(nnti_gni_work_completion));

    wr->wc_post_desc.local_addr            =(uint64_t)&wr->wc;
    wr->wc_post_desc.local_mem_hndl.qword1 =wr->wc_mem_hdl.qword1;
    wr->wc_post_desc.local_mem_hndl.qword2 =wr->wc_mem_hdl.qword2;
    wr->wc_post_desc.remote_addr           =wr->wc_dest_addr;
    wr->wc_post_desc.remote_mem_hndl.qword1=wr->wc_dest_mem_hdl.qword1;
    wr->wc_post_desc.remote_mem_hndl.qword2=wr->wc_dest_mem_hdl.qword2;
    wr->wc_post_desc.length                =sizeof(nnti_gni_work_completion);

    print_raw_buf((void *)wr->wc_post_desc.local_addr, wr->wc_post_desc.length);

    gni_connection *conn=get_conn_peer(peer_hdl);
    assert(conn);
    ep_hdl=conn->ep_hdl;

    GNI_EpSetEventData(
            ep_hdl,
            hash6432shift((uint64_t)src_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf),
            hash6432shift((uint64_t)dest_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.buf));

    if (config.rdma_mode==RDMA_CROSSOVER) {
        if (wr->wc_post_desc.length > config.fma_bte_crossover_size) {
            use_fma=false;
        } else {
            use_fma=true;
        }
    } else if (config.rdma_mode==RDMA_BTE) {
        use_fma=false;
    } else if ((config.rdma_mode==RDMA_FMA) || (config.rdma_mode==RDMA_MIXED)) {
        use_fma=true;
    }

    if (use_fma) {
        log_debug(nnti_debug_level, "calling PostFma(send_buffer_wc ep_hdl(%llu) wc_cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                ep_hdl, transport_global_data.ep_cq_hdl, wr->wc_post_desc.local_addr, wr->wc_post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostFma(ep_hdl, &wr->wc_post_desc);
        trios_stop_timer("PostFma wc", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostFma(fma put) failed: %d", rc);
    } else {
        log_debug(nnti_debug_level, "calling PostRdma(send_buffer_wc ep_hdl(%llu) wc_cq_hdl(%llu), local_addr=%llu, remote_addr=%llu)",
                ep_hdl, transport_global_data.ep_cq_hdl, wr->wc_post_desc.local_addr, wr->wc_post_desc.remote_addr);
        trios_start_timer(call_time);
        rc=GNI_PostRdma(ep_hdl, &wr->wc_post_desc);
        trios_stop_timer("PostRdma wc", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostRdma(rdma put) failed: %d", rc);
    }

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int buffer_send(
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *src_hdl,
        const NNTI_buffer_t *dest_hdl)
{
    int rc=0;
    uint64_t offset=0;
    trios_declare_timer(call_time);

    gni_memory_handle *gni_mem_hdl=NULL;
    gni_work_request  *wr=NULL;

    uint32_t wc_size=sizeof(nnti_gni_work_completion);

    gni_mem_hdl=(gni_memory_handle *)src_hdl->transport_private;
    assert(gni_mem_hdl);
    if (config.use_wr_pool) {
        wr=wr_pool_initiator_pop();
    } else {
        wr=(gni_work_request *)calloc(1, sizeof(gni_work_request));
    }
    assert(wr);

    wr->reg_buf = src_hdl;

    log_debug(nnti_ee_debug_level, "enter");

    if (!config.use_wr_pool) {
        wr->wc_registered  =FALSE;
    }
    wr->wc.ack_received=0;
    wr->wc.inst_id     =transport_global_data.instance;
    wr->wc.op         =GNI_OP_SEND_BUFFER;
    wr->wc.byte_len   =src_hdl->payload_size;
    wr->wc.byte_offset=0;
    wr->wc.src_offset =0;
    wr->wc.dest_offset=0;

    wr->wc_dest_addr          =dest_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.wc_addr;
    wr->wc_dest_mem_hdl.qword1=dest_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword1;
    wr->wc_dest_mem_hdl.qword2=dest_hdl->buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword2;

    nthread_lock(&gni_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "calling send_buffer()");
    trios_start_timer(call_time);
    rc=send_buffer(peer_hdl, src_hdl, dest_hdl, wr);
    trios_stop_timer("send_buffer", call_time);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "send_buffer() failed: %d", rc);

    log_debug(nnti_debug_level, "calling send_buffer_wc()");
    trios_start_timer(call_time);
    rc=send_buffer_wc(peer_hdl, src_hdl, dest_hdl, wr);
    trios_stop_timer("send_buffer_wc", call_time);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "send_buffer_wc() failed: %d", rc);

    gni_mem_hdl->wr_queue.push_back(wr);
    insert_wr_wrhash(wr);

    nthread_unlock(&gni_mem_hdl->wr_queue_lock);

    wr->last_op      =GNI_OP_SEND_BUFFER;
    wr->peer_instance=dest_hdl->buffer_owner.peer.NNTI_remote_process_t_u.gni.inst_id;
    log_debug(nnti_event_debug_level, "wr->last_peer==%lu", wr->peer_instance);

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}

static int send_unblock(
        gni_request_queue_handle *local_req_queue_attrs)
{
    int rc=0;
    gni_post_descriptor_t  post_desc;
    uint32_t *ptr32=NULL;

    trios_declare_timer(call_time);

    log_debug(nnti_ee_debug_level, "enter");


    conn_by_inst_iter_t i;
    for (i=connections_by_instance.begin(); i != connections_by_instance.end(); i++) {

        //NNTI_instance_id key = i->first;
        gni_connection *conn = i->second;

        if (conn==NULL) {
            log_error(nnti_debug_level, "connection in instance map is NULL");
            print_instance_map();
        }

        memset(&post_desc, 0, sizeof(gni_post_descriptor_t));
        post_desc.type           =GNI_POST_CQWRITE;
        post_desc.cq_mode        =GNI_CQMODE_GLOBAL_EVENT|GNI_CQMODE_REMOTE_EVENT;

        set_dlvr_mode(&post_desc);
        set_rdma_mode(&post_desc);

        post_desc.remote_mem_hndl=conn->queue_remote_attrs.client.unblock_mem_hdl;
        ptr32=(uint32_t*)&post_desc.cqwrite_value;
        ptr32++;
        *ptr32=local_req_queue_attrs->total_req_processed;


        rc=GNI_EpBind (local_req_queue_attrs->unblock_ep_hdl, conn->peer_alps_info.local_addr, conn->peer_instance);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpBind(reset_ep_hdl, inst_id=%llu) failed: %d", conn->peer_instance, rc);

        log_debug(nnti_debug_level, "calling PostCqWrite(send_unblock to instance(%llu))", conn->peer_instance);
        trios_start_timer(call_time);
        rc=GNI_PostCqWrite(local_req_queue_attrs->unblock_ep_hdl, &post_desc);
        trios_stop_timer("send_unblock", call_time);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "PostCqWrite(send_unblock, inst_id=%llu) failed: %d", conn->peer_instance, rc);

        post_wait(local_req_queue_attrs->unblock_cq_hdl, 1000, 0);

        rc=GNI_EpUnbind (local_req_queue_attrs->unblock_ep_hdl);
        if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpUnbind(inst_id=%llu) failed: %d", conn->peer_instance, rc);
    }

    log_debug(nnti_ee_debug_level, "exit");

    return(0);
}


static int client_req_queue_init(
        gni_connection *c)
{
    int rc;

    nnti_gni_client_queue  *q              =&c->queue_local_attrs;
    alpsAppGni_t           *server_params  =&c->peer_alps_info;
    uint64_t                server_instance=c->peer_instance;

    log_debug(nnti_debug_level, "enter");

    q->last_offset=0;

    q->req_index=0;
    q->req_index_addr=(uint64_t)&q->req_index;

    log_debug(nnti_debug_level, "client_req_queue->req_index_addr=%llu", (uint64_t)q->req_index_addr);

    q->unblock_buffer=0;
    q->unblock_buffer_addr=(uint64_t)&q->unblock_buffer;

    rc=GNI_CqCreate (
            c->nic_hdl,
            1,
            0,
            GNI_CQ_BLOCKING,
            NULL,
            NULL,
            &q->req_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqCreate() failed: %d", rc);
    rc=GNI_CqCreate (
            c->nic_hdl,
            1,
            0,
            GNI_CQ_BLOCKING,
            NULL,
            NULL,
            &q->req_index_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqCreate() failed: %d", rc);
    rc=GNI_CqCreate (
            c->nic_hdl,
            5,
            0,
            GNI_CQ_BLOCKING,
            NULL,
            NULL,
            &q->unblock_mem_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqCreate() failed: %d", rc);

    rc=GNI_MemRegister (
            c->nic_hdl,
            q->req_index_addr,
            sizeof(uint64_t),
            NULL,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &q->req_index_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemRegister(1) failed: %d", rc);
    rc=GNI_MemRegister (
            c->nic_hdl,
            q->unblock_buffer_addr,
            sizeof(uint64_t),
            q->unblock_mem_cq_hdl,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &q->unblock_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemRegister(1) failed: %d", rc);

    rc=GNI_EpCreate (
            c->nic_hdl,
            q->req_cq_hdl,
            &q->req_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpCreate(req_ep_hdl) failed: %d", rc);
    rc=GNI_EpBind (
            q->req_ep_hdl,
            server_params->local_addr,
            server_instance);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpBind(req_ep_hdl) failed: %d", rc);

    rc=GNI_EpCreate (
            c->nic_hdl,
            q->req_index_cq_hdl,
            &q->req_index_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpCreate(req_index_ep_hdl) failed: %d", rc);
    rc=GNI_EpBind (
            q->req_index_ep_hdl,
            server_params->local_addr,
            server_instance);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpBind(req_index_ep_hdl) failed: %d", rc);

    log_debug(nnti_debug_level, "exit");
}

static int client_req_queue_destroy(
        gni_connection *c)
{
    int rc;
    log_level debug_level = nnti_debug_level;  // nnti_debug_level;

    log_debug(debug_level, "enter");

    nnti_gni_client_queue     *q=&c->queue_local_attrs;

    rc=GNI_MemDeregister (c->nic_hdl, &q->req_index_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "MemDeregister(1) failed: %d", rc);
    rc=GNI_MemDeregister (c->nic_hdl, &q->unblock_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "MemDeregister(1) failed: %d", rc);

    rc=GNI_CqDestroy (q->unblock_mem_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "CqDestroy() failed: %d", rc);

    rc=GNI_EpUnbind (q->req_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "EpUnbind() failed: %d", rc);
    rc=GNI_EpDestroy (q->req_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "EpDestroy() failed: %d", rc);

    rc=GNI_EpUnbind (q->req_index_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "EpUnbind() failed: %d", rc);
    rc=GNI_EpDestroy (q->req_index_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "EpDestroy() failed: %d", rc);

    rc=GNI_CqDestroy (q->req_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "CqDestroy() failed: %d", rc);
    rc=GNI_CqDestroy (q->req_index_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(debug_level, "CqDestroy() failed: %d", rc);

    log_debug(debug_level, "exit");

    return(0);
}

static int server_req_queue_init(
        gni_request_queue_handle *q,
        char                     *buffer,
        uint64_t                  req_size,
        uint64_t                  req_count)
{
    int rc;
    gni_memory_handle *gni_mem_hdl=(gni_memory_handle *)q->reg_buf->transport_private;

    log_debug(nnti_debug_level, "enter");

    q->req_buffer     =buffer;
    q->req_size       =req_size;
    q->req_count      =req_count;
    q->req_buffer_size=req_count*req_size;

    q->wc_buffer_size=q->req_count * sizeof(nnti_gni_work_completion);
    q->wc_buffer     =(nnti_gni_work_completion *)calloc(q->req_count, sizeof(nnti_gni_work_completion));

    q->req_index     =0;
    q->req_index_addr=(uint64_t)&q->req_index;

    log_debug(nnti_debug_level, "server_req_queue->req_index_addr=%llu", (uint64_t)q->req_index_addr);

    q->req_processed=0;
    q->req_processed_reset_limit=q->req_count;

    rc=GNI_CqCreate (
            transport_global_data.nic_hdl,
            1,
            0,
            GNI_CQ_BLOCKING,
            NULL,
            NULL,
            &q->unblock_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqCreate() failed: %d", rc);
    rc=GNI_EpCreate (
            transport_global_data.nic_hdl,
            q->unblock_cq_hdl,
            &q->unblock_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpCreate(reset_ep_hdl) failed: %d", rc);

    rc=GNI_CqCreate (
            transport_global_data.nic_hdl,
            req_count,
            0,
            GNI_CQ_BLOCKING,
            NULL,
            NULL,
            &q->wc_mem_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqCreate() failed: %d", rc);

    rc=GNI_MemRegister (
            transport_global_data.nic_hdl,
            (uint64_t)q->req_buffer,
            q->req_buffer_size,
            NULL,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &gni_mem_hdl->mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemRegister(1) failed: %d", rc);
    rc=GNI_MemRegister (
            transport_global_data.nic_hdl,
            q->req_index_addr,
            sizeof(uint64_t),
            NULL,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &q->req_index_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemRegister(1) failed: %d", rc);
    rc=GNI_MemRegister (
            transport_global_data.nic_hdl,
            (uint64_t)q->wc_buffer,
            q->wc_buffer_size,
            q->wc_mem_cq_hdl,
            GNI_MEM_READWRITE,
            (uint32_t)-1,
            &q->wc_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemRegister(1) failed: %d", rc);

    log_debug(nnti_debug_level, "exit");
}

static int server_req_queue_destroy(
        gni_request_queue_handle *q)
{
    int rc;
    gni_memory_handle *gni_mem_hdl=(gni_memory_handle *)q->reg_buf->transport_private;

    log_debug(nnti_debug_level, "enter");

    rc=GNI_MemDeregister (
            transport_global_data.nic_hdl,
            &gni_mem_hdl->mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemDeregister(1) failed: %d", rc);
    rc=GNI_MemDeregister (
            transport_global_data.nic_hdl,
            &q->req_index_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemDeregister(1) failed: %d", rc);
    rc=GNI_MemDeregister (
            transport_global_data.nic_hdl,
            &q->wc_mem_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "MemDeregister(1) failed: %d", rc);

    rc=GNI_CqDestroy (q->wc_mem_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqDestroy() failed: %d", rc);

    rc=GNI_EpDestroy (q->unblock_ep_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "EpDestroy() failed: %d", rc);
    rc=GNI_CqDestroy (q->unblock_cq_hdl);
    if (rc!=GNI_RC_SUCCESS) log_error(nnti_debug_level, "CqDestroy() failed: %d", rc);

    free(q->wc_buffer);

    log_debug(nnti_debug_level, "exit");
}

#define LISTEN_IFACE_BASENAME "ipog"
static NNTI_result_t get_ipaddr(
        char *ipaddr,
        int maxlen)
{
    struct ifaddrs * ifAddrStruct=NULL;
    struct ifaddrs * ifa=NULL;
    void * tmpAddrPtr=NULL;

    getifaddrs(&ifAddrStruct);
    for (ifa = ifAddrStruct; ifa != NULL; ifa = ifa->ifa_next) {
        ipaddr[0]='\0';
        if (ifa ->ifa_addr->sa_family==AF_INET) { // check it is IP4
            // is a valid IP4 Address
            tmpAddrPtr=&((struct sockaddr_in *)ifa->ifa_addr)->sin_addr;
            log_debug(nnti_debug_level, "checking iface (IPv4) name (%s)", ifa->ifa_name);
            inet_ntop(AF_INET, tmpAddrPtr, ipaddr, maxlen);
            log_debug(nnti_debug_level, "hostname(%s) has IP Address %s", ifa->ifa_name, ipaddr);
            if (0==strncmp(ifa->ifa_name, LISTEN_IFACE_BASENAME, strlen(LISTEN_IFACE_BASENAME))) {
                log_debug(nnti_debug_level, "hostname(%s) matches", ifa->ifa_name);
                break;
            }
        } else if (ifa->ifa_addr->sa_family==AF_INET6) { // check it is IP6
            // is a valid IP6 Address
            tmpAddrPtr=&((struct sockaddr_in *)ifa->ifa_addr)->sin_addr;
            log_debug(nnti_debug_level, "checking iface (IPv6) name (%s)", ifa->ifa_name);
            inet_ntop(AF_INET6, tmpAddrPtr, ipaddr, maxlen);
            log_debug(nnti_debug_level, "hostname(%s) has IP Address %s", ifa->ifa_name, ipaddr);
            if (0==strncmp(ifa->ifa_name, LISTEN_IFACE_BASENAME, strlen(LISTEN_IFACE_BASENAME))) {
                log_debug(nnti_debug_level, "hostname(%s) matches", ifa->ifa_name);
                break;
            }
        }
    }
    if (ifAddrStruct!=NULL) freeifaddrs(ifAddrStruct);

    log_debug(nnti_debug_level, "ipaddr(%s)", ipaddr);
    if (ipaddr[0]=='\0') {
        return NNTI_ENOENT;
    } else {
        return NNTI_OK;
    }
}

static void config_init(nnti_gni_config *c)
{
    c->use_alps_ptag                    =true;
    c->use_wr_pool                      =false;
    c->use_rdma_target_ack              =true;
    c->use_rdma_events                  =true;
    c->use_rdma_fence                   =false;
    c->pi_ordering                      =PI_ORDERING_DEFAULT; /* DEFAULT, STRICT, RELAXED */
    c->rdma_mode                        =RDMA_MIXED;          /* FMA, BTE, MIXED, CROSSOVER */
    c->fma_bte_crossover_size           =4096;
}

static void config_get_from_env(nnti_gni_config *c)
{
    char *env_str=NULL;

    if ((env_str=getenv("TRIOS_NNTI_USE_ALPS_PTAG")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->use_alps_ptag to TRUE");
            c->use_alps_ptag=true;
        } else {
            log_debug(nnti_debug_level, "setting c->use_alps_ptag to FALSE");
            c->use_alps_ptag=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_ALPS_PTAG is undefined.  using c->use_alps_ptag default");
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
    if ((env_str=getenv("TRIOS_NNTI_USE_RDMA_EVENTS")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->use_rdma_events to TRUE");
            c->use_rdma_events=true;
        } else {
            log_debug(nnti_debug_level, "setting c->use_rdma_events to FALSE");
            c->use_rdma_events=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_RDMA_EVENTS is undefined.  using c->use_rdma_events default");
    }
    if ((env_str=getenv("TRIOS_NNTI_USE_RDMA_FENCE")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->use_rdma_fence to TRUE");
            c->use_rdma_fence=true;
        } else {
            log_debug(nnti_debug_level, "setting c->use_rdma_fence to FALSE");
            c->use_rdma_fence=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_RDMA_FENCE is undefined.  using c->use_rdma_fence default");
    }
    if ((env_str=getenv("TRIOS_NNTI_PI_ORDERING")) != NULL) {
        if (!strcasecmp(env_str, "PI_ORDERING_DEFAULT")) {
            log_debug(nnti_debug_level, "setting c->pi_ordering to PI_ORDERING_DEFAULT");
            c->pi_ordering=PI_ORDERING_DEFAULT;
        } else if (!strcasecmp(env_str, "PI_ORDERING_STRICT")) {
            log_debug(nnti_debug_level, "setting c->pi_ordering to PI_ORDERING_STRICT");
            c->pi_ordering=PI_ORDERING_STRICT;
        } else if (!strcasecmp(env_str, "PI_ORDERING_RELAXED")) {
            log_debug(nnti_debug_level, "setting c->pi_ordering to PI_ORDERING_RELAXED");
            c->pi_ordering=PI_ORDERING_RELAXED;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_PI_ORDERING is undefined.  using c->pi_ordering default");
    }
    if ((env_str=getenv("TRIOS_NNTI_RDMA_MODE")) != NULL) {
        if (!strcasecmp(env_str, "RDMA_FMA")) {
            log_debug(nnti_debug_level, "setting c->rdma_mode to RDMA_FMA");
            c->rdma_mode=RDMA_FMA;
        } else if (!strcasecmp(env_str, "RDMA_BTE")) {
            log_debug(nnti_debug_level, "setting c->rdma_mode to RDMA_BTE");
            c->rdma_mode=RDMA_BTE;
        } else if (!strcasecmp(env_str, "RDMA_MIXED")) {
            log_debug(nnti_debug_level, "setting c->rdma_mode to RDMA_MIXED");
            c->rdma_mode=RDMA_MIXED;
        } else if (!strcasecmp(env_str, "RDMA_CROSSOVER")) {
            log_debug(nnti_debug_level, "setting c->rdma_mode to RDMA_CROSSOVER");
            c->rdma_mode=RDMA_CROSSOVER;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_RDMA_MODE is undefined.  using c->rdma_mode default");
    }
    if ((env_str=getenv("TRIOS_NNTI_FMA_BTE_CROSSOVER_SIZE")) != NULL) {
        errno=0;
        uint32_t crossover_size=strtoul(env_str, NULL, 0);
        if (errno == 0) {
            log_debug(nnti_debug_level, "setting c->fma_bte_crossover_size to %lu", crossover_size);
            c->fma_bte_crossover_size=crossover_size;
        } else {
            log_debug(nnti_debug_level, "TRIOS_NNTI_FMA_BTE_CROSSOVER_SIZE value conversion failed (%s).  using c->fma_bte_crossover_size default.", strerror(errno));
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_FMA_BTE_CROSSOVER_SIZE is undefined.  using c->fma_bte_crossover_size default");
    }
}
