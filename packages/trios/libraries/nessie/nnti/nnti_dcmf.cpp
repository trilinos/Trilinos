/*
 * nnti_dcmf.c
 *
 *  Created on: April 19, 2012
 */

#include "Trios_config.h"
#include "Trios_threads.h"
#include "Trios_timer.h"
#include "Trios_signal.h"
#include "Trios_nnti_fprint_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <spi/bgp_SPI.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>

#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include <sys/time.h>
#include <sys/resource.h>

#include <map>
#include <deque>

#include <dcmf.h>
#include "nnti_dcmf.h"
#include "nnti_utils.h"
#include "nnti_internal.h"

#define TRACE_ERR(x)		/* fprintf  x */

/* GLOBAL VAR */

DCMF_Protocol_t send_prot;
DCMF_Protocol_t control0_prot;
DCMF_Control_t info;
DCMF_Protocol_t req_index_prot;

DCMF_Send_Configuration_t conf;	/* Register the send protocol for intial conn */

#define BUFSIZE  2048
#define ITERATIONS  1000
#define MAX_CONNECTION  1024
size_t my_rank;
volatile unsigned _recv_active;
volatile unsigned send_active;
volatile unsigned _recv_iteration;
volatile unsigned confirm_conn;
DCMF_Request_t _recv_request[ITERATIONS];
char _recv_buffer[BUFSIZE] __attribute__ ((__aligned__ (16)));
volatile static int _is_server;

/*** End of GLOBAL VAR ***/


/**
 * These states are used to signal events between the completion handler
 * and the main client or server thread.
 *
 * Once CONNECTED, they cycle through RDMA_READ_ADV, RDMA_WRITE_ADV,
 * and RDMA_WRITE_COMPLETE for each ping.
 */
typedef enum
{
  IDLE = 1,
  CONNECT_REQUEST,
  CONNECTED,
  DISCONNECTED,
  ERROR
} bgpdcmf_connection_state;

typedef enum
{
  SERVER_CONNECTION,
  CLIENT_CONNECTION
} bgpdcmf_connection_type;

typedef enum
{
  REQUEST_BUFFER,
  RESULT_BUFFER,
  RECEIVE_BUFFER,
  SEND_BUFFER,
  GET_SRC_BUFFER,
  GET_DST_BUFFER,
  PUT_SRC_BUFFER,
  PUT_DST_BUFFER,
  RDMA_TARGET_BUFFER,
  UNKNOWN_BUFFER
} bgpdcmf_buffer_type;


#define DCMF_OP_PUT_INITIATOR  1
#define DCMF_OP_GET_INITIATOR  2
#define DCMF_OP_PUT_TARGET     3
#define DCMF_OP_GET_TARGET     4
#define DCMF_OP_SEND_REQ           5
#define DCMF_OP_NEW_REQUEST    6
#define DCMF_OP_RESULT         7
#define DCMF_OP_RECEIVE        8
#define DCMF_OP_REGISTER_RDMA        9

typedef enum
{
  SEND_COMPLETE,
  RECV_COMPLETE,
  RDMA_WRITE_INIT,
  RDMA_TARGET_INIT,
  RDMA_READ_INIT,
  RDMA_COMPLETE
} bgpdcmf_op_state_t;


typedef struct
{

  size_t req_rank;
  DCMF_Memregion_t response_hdl;
} nnti_bgpdcmf_fetch_index_req ALIGN_QUADWORD;

typedef struct
{
  uint32_t req_received;
  int local_send_complete;
  int local_get_complete;
  int remote_put_complete;
  int op_complete;
  uint32_t req_rank;
  uint64_t byte_len;
  uint64_t byte_offset;
  uint64_t src_offset;
  uint64_t dest_offset;
  uint16_t op;
} nnti_bgpdcmf_work_completion ALIGN_QUADWORD;

/**
 * attrs to send to the client
 */
typedef struct
{
  uint64_t req_index_addr;	/* address of the request index var on the server */
  DCMF_Memregion_t req_index_mem_hdl;

  uint64_t req_buffer_addr;	/* address of the request buffer */
  uint64_t req_size;		/* the maximum size of a request */
  uint64_t req_count;		/* the number of requests that will fit in the queue */
  DCMF_Memregion_t req_mem_hdl;
  uint64_t wc_buffer_addr;	/* address of the work completion array */
  DCMF_Memregion_t wc_mem_hdl;
  uint64_t unblock_buffer_addr;	/* address of the fetch_req_array  */
  DCMF_Memregion_t unblock_mem_hdl;
} nnti_bgpdcmf_server_queue_attrs;

/**
 * attrs to send to the server
 */
typedef struct
{
  uint64_t req_index;		/* this buffer contains the current index of the request queue */
  uint64_t req_index_addr;	/* address of the index buffer */
  DCMF_Memregion_t req_index_mem_hdl;	/* mem hdl to send to the server */
  DCMF_Memregion_t fetch_index_hdl;
} nnti_bgpdcmf_client_queue_attrs;

typedef union
{
  nnti_bgpdcmf_client_queue_attrs client;
  nnti_bgpdcmf_server_queue_attrs server;
} nnti_bgpdcmf_queue_remote_attrs;

typedef struct
{
  uint32_t peer_rank;
  nnti_bgpdcmf_client_queue_attrs queue_local_attrs;
  nnti_bgpdcmf_queue_remote_attrs queue_remote_attrs;
  bgpdcmf_connection_state state;
  bgpdcmf_connection_type connection_type;
} bgpdcmf_connection;


typedef struct
{
  uint8_t is_initiator;

  const NNTI_buffer_t *reg_buf;
  nnti_bgpdcmf_work_completion wc;
  uint64_t wc_dest_addr;
  DCMF_Memregion_t wc_dest_mem_hdl;
  DCMF_Memregion_t wc_mem_hdl;
  uint8_t last_op;
  bgpdcmf_op_state_t op_state;
  uint64_t is_last_op_complete;
  uint32_t peer_rank;
} dcmf_work_request;


typedef
  std::deque <
dcmf_work_request * >
  wr_queue_t;
typedef
  std::deque <
dcmf_work_request * >::iterator
  wr_queue_iter_t;

typedef struct
{
  bgpdcmf_buffer_type
    type;
  DCMF_Memregion_t
    mem_hdl;			/* actual data memory handle */
  wr_queue_t
    wr_queue;
  nthread_lock_t
    wr_queue_lock;
  uint32_t
    ref_count;
} bgpdcmf_memory_handle;

typedef struct
{
  NNTI_buffer_t *
    reg_buf;

  uint64_t
    last_index_before_reset;
  uint64_t
    req_index;			/* index of the next available slot in the request queue */
  uint64_t
    req_index_addr;		/* address of the request index var on the server */
  DCMF_Memregion_t
    req_index_mem_hdl;

  char *
    req_buffer;			/* pointer to the head of the request buffer */
  uint64_t
    req_size;			/* the maximum size of a request */
  uint64_t
    req_count;			/* the number of requests that will fit in the queue */
  uint64_t
    req_buffer_size;		/* the size of the request buffer in bytes (req_size*req_count) */
  DCMF_Memregion_t
    req_mem_hdl;

  nnti_bgpdcmf_work_completion *
    wc_buffer;			/* pointer to work completion array */
  uint64_t
    wc_buffer_size;		/* the size of the work completion buffer in bytes  */
  DCMF_Memregion_t
    wc_mem_hdl;
  nnti_bgpdcmf_fetch_index_req *
    unblock_buffer;		/* pointer to individual fetch index array */
  DCMF_Memregion_t
    unblock_mem_hdl;
  uint64_t
    req_processed_reset_limit;
  uint64_t
    req_processed;
  uint64_t
    total_req_processed;
} bgpdcmf_request_queue_handle;

typedef struct
{
  uint32_t
    myrank;
  uint32_t
    remote_rank;
  int
    mypid;
  size_t  conn_req[1024];
  bgpdcmf_request_queue_handle
    req_queue;
} bgpdcmf_transport_global;

static nthread_lock_t
  nnti_bgpdcmf_lock;
static nthread_lock_t
  nnti_index_lock;
static int
start_connection_listener_thread (void);
static int
start_index_thread (void);
static NNTI_result_t
register_memory (bgpdcmf_memory_handle * hdl,
		 void *buf, uint64_t len, const NNTI_buf_ops_t remote_op);
static int
unregister_memory (bgpdcmf_memory_handle * hdl);
static NNTI_result_t
process_event (const NNTI_buffer_t * reg_buf,
	       const NNTI_buf_ops_t remote_op,
	       dcmf_work_request * wr, const int timeout);
/* NEW SHYAMALI */

static NNTI_result_t
repost_recv_work_request (NNTI_buffer_t * reg_buf, dcmf_work_request * wr);

static dcmf_work_request *
first_incomplete_wr (bgpdcmf_memory_handle * _hdl);
static int8_t
is_buf_op_complete (const NNTI_buffer_t * reg_buf);
static int8_t
is_wr_complete (dcmf_work_request * wr);
static int8_t
is_any_buf_op_complete (const NNTI_buffer_t ** buf_list,
			const uint32_t buf_count, uint32_t * which);
static int8_t
is_all_buf_ops_complete (const NNTI_buffer_t ** buf_list,
			 const uint32_t buf_count);
static void
create_status (const NNTI_buffer_t * reg_buf,
	       const NNTI_buf_ops_t remote_op,
	       int nnti_rc, NNTI_status_t * status);

/* NEW */
static void
create_peer (NNTI_peer_t * peer, int rank);
static NNTI_result_t
inject_full (int rank, const void *buf, size_t num);
static NNTI_result_t
init_connection (bgpdcmf_connection ** conn, int rank, const int is_server);
static void
close_connection (bgpdcmf_connection * c);

static NNTI_result_t
insert_conn_rank (const uint32_t rank, bgpdcmf_connection * conn);
static bgpdcmf_connection *
get_conn_rank (const uint32_t rank);
static bgpdcmf_connection *
del_conn_rank (const uint32_t rank);

static int
server_req_queue_init (bgpdcmf_request_queue_handle * q,
		       char *buffer, uint64_t req_size, uint64_t req_count);
static int
server_req_queue_destroy (bgpdcmf_request_queue_handle * q);

static int
client_req_queue_init (bgpdcmf_connection * c);
static int
client_req_queue_destroy (bgpdcmf_connection * c);

static int
reset_req_index (bgpdcmf_request_queue_handle * req_queue_attrs);

static int
fetch_server_req_buffer_offset (nnti_bgpdcmf_client_queue_attrs *
				local_req_queue_attrs,
				nnti_bgpdcmf_server_queue_attrs *
				remote_req_queue_attrs,
				uint64_t addend, uint64_t * offset, int rank);
static int
send_req (nnti_bgpdcmf_client_queue_attrs * local_req_queue_attrs,
	  nnti_bgpdcmf_server_queue_attrs * remote_req_queue_attrs,
	  uint64_t offset, const NNTI_buffer_t * reg_buf,
	  int rank, dcmf_work_request * wr);

static NNTI_result_t
send_req_wc (nnti_bgpdcmf_client_queue_attrs * local_req_queue_attrs,
	     nnti_bgpdcmf_server_queue_attrs * remote_req_queue_attrs,
	     uint64_t offset,
	     const NNTI_buffer_t * reg_buf, dcmf_work_request * wr, int rank);

static int
request_send (nnti_bgpdcmf_client_queue_attrs * client_q,
	      nnti_bgpdcmf_server_queue_attrs * server_q,
	      const NNTI_buffer_t * reg_buf, int rank);

/* NEW */
static NNTI_result_t
register_wc (dcmf_work_request * wr);
static NNTI_result_t
insert_buf_bufhash (NNTI_buffer_t * buf);
static NNTI_buffer_t *
get_buf_bufhash (const uint32_t bufhash);
static NNTI_buffer_t *
del_buf_bufhash (NNTI_buffer_t * buf);
static void
print_bufhash_map (void);

static NNTI_result_t
insert_wr_wrhash (dcmf_work_request *);
static dcmf_work_request *
get_wr_wrhash (const uint32_t bufhash);
static dcmf_work_request *
del_wr_wrhash (dcmf_work_request *);
static void
print_wrhash_map (void);
static void
send_rdma_wc (dcmf_work_request * wr,
	      const NNTI_buffer_t * local_buf,
	      const NNTI_buffer_t * remote_buf);



/* Thomas Wang's 64 bit to 32 bit Hash Function (http://www.concentric.net/~ttwang/tech/inthash.htm) */
static uint32_t
hash6432shift (uint64_t key)
{
  key = (~key) + (key << 18);	// key = (key << 18) - key - 1;
  key = key ^ (key >> 31);
  key = key * 21;		// key = (key + (key << 2)) + (key << 4);
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);
  return (uint32_t) key;
}

/*
 * We need a couple of maps to keep track of connections.  Servers need to find
 * connections by QP number when requests arrive.  Clients need to find connections
 * by peer address and port.  Setup those maps here.
 */
static
  std::map < int,
bgpdcmf_connection * >connections_by_rank;
typedef
  std::map < int,
bgpdcmf_connection * >::iterator
  conn_by_rank_iter_t;
typedef
  std::pair < int,
bgpdcmf_connection * >
  conn_by_rank_t;
static nthread_lock_t
  nnti_conn_peer_lock;

static int8_t
is_wr_queue_empty (const NNTI_buffer_t * buf);

static
  std::map <
uint32_t, NNTI_buffer_t * >buffers_by_bufhash;
typedef
  std::map <
  uint32_t,
NNTI_buffer_t * >::iterator
  buf_by_bufhash_iter_t;
typedef
  std::pair <
  uint32_t,
NNTI_buffer_t * >
  buf_by_bufhash_t;
static nthread_lock_t
  nnti_buf_bufhash_lock;

static
  std::map <
  uint32_t,
dcmf_work_request * >
  wr_by_wrhash;
typedef
  std::map <
  uint32_t,
dcmf_work_request * >::iterator
  wr_by_wrhash_iter_t;
typedef
  std::pair <
  uint32_t,
dcmf_work_request * >
  wr_by_wrhash_t;
static nthread_lock_t
  nnti_wr_wrhash_lock;



/* NEW */

static bgpdcmf_transport_global
  transport_global_data;
static const int
  MIN_TIMEOUT = 100000000;		/* in milliseconds */
static log_level
  nnti_event_debug_level;
static nnti_bgpdcmf_fetch_index_req
  client_fetch_req;

/**************************  Callback functions for DCMF_Send / Recv protocol ************/
/* Low overhead timing function */

typedef unsigned long long
  ticks;

static __inline__ ticks
getticks (void)
{
  unsigned int
    tbl,
    tbu0,
    tbu1;

  do
    {
      __asm__
      __volatile__ ("mftbu %0":"=r" (tbu0));
      __asm__
      __volatile__ ("mftb %0":"=r" (tbl));
      __asm__
      __volatile__ ("mftbu %0":"=r" (tbu1));
    }
  while (tbu0 != tbu1);

  return (((unsigned long long) tbu0) << 32) | tbl;
}


void *
allocateMemory (size_t size)
{
  void *
    ptr;
  int
    rc;

  rc = posix_memalign (&ptr, 32, size);
  if (rc)
    return NULL;

  return ptr;
}


static void
decrement (void *clientdata, DCMF_Error_t * err)
{
  unsigned *
    _clientdata = (unsigned *) clientdata;
  --*_clientdata;
  TRACE_ERR ((stderr, "(%zd) decrement() clientdata = %p, %d => %d\n",
	      my_rank, _clientdata, *_clientdata, *_clientdata - 1));
}


void control_recv (void                 * clientdata,
                   const DCMF_Control_t * info,
                   size_t                 peer)
{
  DCMF_Control_t tmp;
  memcpy((void *)&tmp, (const void *)info, sizeof(DCMF_Control_t));
  (*((unsigned *)clientdata))--;
  //log_debug(nnti_debug_level, "Received a confirmation on connection\n");
}

/* --------------------------------------------------------------- */

static void
cb_recv_new_short (void *clientdata,
		   const DCQuad * msginfo,
		   unsigned count, size_t peer, const char *src, size_t bytes)
{
  unsigned *
    _clientdata = (unsigned *) clientdata;
  NNTI_result_t
    rc = NNTI_OK;

  bgpdcmf_connection *
    conn = NULL;

  memcpy (&_recv_buffer, src, bytes);
  nthread_lock (&nnti_bgpdcmf_lock);
  if (strncmp (_recv_buffer, "INIT", 4) == 0)
    {
       conn = get_conn_rank(peer);
       if (conn == NULL) {
		transport_global_data.conn_req[peer] = peer;
       }
    }
    else if(_is_server == 1){
	struct
  	{
    		nnti_bgpdcmf_client_queue_attrs
      		client_attrs;
  	} ca_in;
  	memset (&ca_in, 0, sizeof (ca_in));
        memcpy(&ca_in, src, sizeof (ca_in));
	conn = get_conn_rank(peer);
	 if (conn != NULL)
		conn->queue_remote_attrs.client = ca_in.client_attrs;
	  DCMF_CriticalSection_enter (0);
          DCMF_Control (&control0_prot, DCMF_MATCH_CONSISTENCY, peer, &info);
	  DCMF_CriticalSection_exit (0);
  // log_debug(nnti_debug_level, "(%zd) cb_recv_new_short() recv full client attr  bytes  %d : clientdata = %p, %d => %d\n", my_rank, bytes, _clientdata, *_clientdata, *_clientdata-1);
    }
  --*_clientdata;
    nthread_unlock (&nnti_bgpdcmf_lock);
}

/* -------------------------------------------------------------- */

static DCMF_Request_t *
cb_recv_new (void *clientdata,
	     const DCQuad * msginfo,
	     unsigned count,
	     size_t senderID,
	     size_t sndlen,
	     size_t * rcvlen, char **rcvbuf, DCMF_Callback_t * cb_info)
{
  TRACE_ERR ((stderr, "(%zd) cb_recv_new() clientdata = %p, sndlen = %zd\n",
	      my_rank, clientdata, sndlen));
  *rcvbuf = (char *) _recv_buffer;
  *rcvlen = sndlen > 0 ? sndlen : 1;
  cb_info->function = decrement;
  cb_info->clientdata = clientdata;

  // log_debug(nnti_debug_level, "cb_recv_new long case\n");
  return &_recv_request[_recv_iteration++];
}



void
send_once (char *buffer, size_t sndlen, size_t targetrank,
	   DCMF_Consistency consistency, DCMF_Protocol_t * protocol)
{
  DCQuad
    msginfo;
  DCMF_Request_t
    request;

  send_active = 1;
  DCMF_Callback_t
  cb_info = { decrement, (void *) &send_active };
  nthread_lock (&nnti_bgpdcmf_lock);
  DCMF_CriticalSection_enter (0);
  DCMF_Send (protocol,
	     &request,
	     cb_info, consistency, targetrank, sndlen, buffer, &msginfo, 1);

  while (send_active)
    DCMF_Messager_advance ();
  DCMF_CriticalSection_exit (0);
  nthread_unlock (&nnti_bgpdcmf_lock);
  //log_debug(nnti_debug_level, "Send once completed %d\n", transport_global_data.myrank);
  send_active = 1;
}

int
dump_memregion (DCMF_Memregion_t * memreg)
{

  fprintf (stderr, "memregion=0x%08x %08x %08x %08x \n",
	   *((unsigned *) (((char *) (memreg)) + 0)),
	   *((unsigned *) (((char *) (memreg)) + 4)),
	   *((unsigned *) (((char *) (memreg)) + 8)),
	   *((unsigned *) (((char *) (memreg)) + 12)));
  return 0;
}


/* ==================== */



/**
 * @brief Initialize NNTI to use a specific transport.
 *
 */
NNTI_result_t
NNTI_bgpdcmf_init (const NNTI_transport_id_t trans_id,
		   const char *my_url, NNTI_transport_t * trans_hdl)
{
  static int
    initialized = 0;

  NNTI_result_t
    rc = NNTI_OK;
  int i;

#ifndef HAS_MPI
  DCMF_Messager_initialize ();
#endif

  my_rank = DCMF_Messager_rank ();

  assert (trans_hdl);

  initialized = 0;
  _recv_active = 1;
  send_active = 1;
  confirm_conn = 1;
  _recv_iteration = 0;

  if (!initialized)
    {
      conf.protocol = DCMF_EAGER_SEND_PROTOCOL;
      conf.network = DCMF_TORUS_NETWORK;
      conf.cb_recv_short = cb_recv_new_short;
      conf.cb_recv_short_clientdata = (void *) &_recv_active;
      conf.cb_recv = cb_recv_new;
      conf.cb_recv_clientdata = (void *) &_recv_active;
      DCMF_Result
	result = DCMF_Send_register (&send_prot, &conf);
      if (result != DCMF_SUCCESS)
	{
	  fprintf (stderr, "DCMF_SEND register failed in _init\n");
	}

      /*   Register DCMF control */
      DCMF_Control_Configuration_t c0_conf;
        c0_conf.protocol = DCMF_DEFAULT_CONTROL_PROTOCOL;
        c0_conf.network = DCMF_TORUS_NETWORK;
  	c0_conf.cb_recv = control_recv;
  	c0_conf.cb_recv_clientdata = (void *)&confirm_conn; 
      result =  DCMF_Control_register(&control0_prot, &c0_conf);

      if (result != DCMF_SUCCESS)
        {
          fprintf (stderr, "DCMF_Control register failed register failed in _init\n");
        }

      memset (&transport_global_data, 0, sizeof (bgpdcmf_transport_global));

      nthread_lock_init (&nnti_bgpdcmf_lock);
      nthread_lock_init (&nnti_index_lock);
      nthread_lock_init (&nnti_conn_peer_lock);
      nthread_lock_init (&nnti_wr_wrhash_lock);
      nthread_lock_init (&nnti_buf_bufhash_lock);

      //log_debug (nnti_debug_level, "my_url=%s", my_url);


      //log_debug (nnti_debug_level, "initializing Blue Gene  DMA DCMF device");

      transport_global_data.mypid = getpid ();
      transport_global_data.myrank = my_rank;
      for (i = 0; i< 1024; ++i)
        {
                transport_global_data.conn_req[i] = -1;
        }

      create_peer (&trans_hdl->me, transport_global_data.myrank);

      initialized = TRUE;
    }

  //log_debug (nnti_debug_level, "exit");

  return (rc);
}


/**
 * @brief Return the URL field of this transport.
 *
 * Return the URL field of this transport.  After initialization, the transport will
 * have a specific location on the network where peers can contact it.  The
 * transport will convert this location to a string that other instances of the
 * transport will recobgpdcmfze.
 *
 */
NNTI_result_t
NNTI_bgpdcmf_get_url (const NNTI_transport_t * trans_hdl,
		      char *url, const uint64_t maxlen)
{
  NNTI_result_t
    rc = NNTI_OK;

  assert (trans_hdl);
  assert (url);
  assert (maxlen > 0);

  //log_debug (nnti_debug_level, "enter");

  strncpy (url, trans_hdl->me.url, maxlen);
  url[maxlen - 1] = '\0';

  //log_debug (nnti_debug_level, "exit");

  return (rc);
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
NNTI_result_t
NNTI_bgpdcmf_connect (const NNTI_transport_t * trans_hdl,
		      const char *url,
		      const int timeout, NNTI_peer_t * peer_hdl)
{
  NNTI_result_t
    rc = NNTI_OK;
  char *
    sendbuf;

  char
    transport[NNTI_URL_LEN];
  char
    address[NNTI_URL_LEN];
  char *
    tmp;
  char *
    temp;
  int
    pset_rank;
  size_t
    sndlen;

  bgpdcmf_connection *
    conn = NULL;


  assert (trans_hdl);
  assert (peer_hdl);

  //log_debug (nnti_debug_level, "enter URL = %s  myrank = %d", url, transport_global_data.myrank);

  if (url != NULL)
    {
      if ((rc =
	   nnti_url_get_transport (url, transport, NNTI_URL_LEN)) != NNTI_OK)
	{
	  return (rc);
	}
      if (0 != strcmp (transport, "dcmf"))
	{
	  /* the peer described by 'url' is not a BGP Torus peer */
	  return (NNTI_EINVAL);
	}

      if ((rc = nnti_url_get_address (url, address, NNTI_URL_LEN)) != NNTI_OK)
	{
	  return (rc);
	}
      pset_rank = strtol (address , NULL, 10);
      transport_global_data.remote_rank = pset_rank;	/* This is server rank */
      //log_debug (nnti_debug_level, "connect client  to %d  ", pset_rank);
    }
  else
    {
      /*  */
      rc = NNTI_EINVAL;
      goto cleanup;
    }
  sendbuf = (char *) allocateMemory (256);	/* send INIT_CONNECTION wth rank  as SW_ARG */
  memset (sendbuf, 0x00, sizeof (sendbuf));
  sprintf (sendbuf, "INIT_CONNECTION:%d:", transport_global_data.myrank);
  sndlen = 32;
  send_once (sendbuf, sndlen, transport_global_data.remote_rank,
	     DCMF_SEQUENTIAL_CONSISTENCY, &send_prot);
  rc = init_connection (&conn, pset_rank, 0);
  if (conn == NULL)
    {
      fprintf (stderr, "bgpdcmf connect error\n");
      rc = NNTI_EIO;
      goto cleanup;
    }
  create_peer (peer_hdl, pset_rank);
  insert_conn_rank (pset_rank, conn);
cleanup:
  free (sendbuf);
  //log_debug (nnti_debug_level, "exit");
  return (rc);
}


/**
 * @brief Terminate communication with this peer.
 *
 * Perform any transport specific actions necessary to end communication with
 * this peer.
 */
NNTI_result_t
NNTI_bgpdcmf_disconnect (const NNTI_transport_t * trans_hdl,
			 NNTI_peer_t * peer_hdl)
{
  NNTI_result_t
    rc = NNTI_OK;

  assert (trans_hdl);
  assert (peer_hdl);

  //log_debug (nnti_debug_level, "enter");
  bgpdcmf_connection *
    conn =
    get_conn_rank (peer_hdl->peer.NNTI_remote_process_t_u.bgpdcmf.pset_rank);
  close_connection (conn);
  del_conn_rank (peer_hdl->peer.NNTI_remote_process_t_u.bgpdcmf.pset_rank);

  free (peer_hdl->url);

  //log_debug (nnti_debug_level, "exit");

  return (rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 */
NNTI_result_t
NNTI_bgpdcmf_register_memory (const NNTI_transport_t * trans_hdl,
			      char *buffer,
			      const uint64_t size,
			      const uint64_t num_elements,
			      const NNTI_buf_ops_t ops,
			      const NNTI_peer_t * peer,
			      NNTI_buffer_t * reg_buf)
{
  NNTI_result_t
    rc = NNTI_OK;

  dcmf_work_request *
    wr = NULL;
  NNTI_buffer_t *
    old_buf = NULL;
  bgpdcmf_memory_handle *
    bgpdcmf_mem_hdl = NULL;


  assert (trans_hdl);
  assert (buffer);
  assert (size > 0);
  assert (ops > 0);
  assert (reg_buf);

  log_debug (nnti_debug_level, "register memory enter buffer=%p", reg_buf);

  old_buf = get_buf_bufhash (hash6432shift ((uint64_t) buffer));
  if (old_buf == NULL)
    {
      bgpdcmf_mem_hdl = new bgpdcmf_memory_handle ();
      bgpdcmf_mem_hdl->ref_count = 1;
    }
  else
    {
      bgpdcmf_mem_hdl = (bgpdcmf_memory_handle *) old_buf->transport_private;
      bgpdcmf_mem_hdl->ref_count++;
    }

  assert (bgpdcmf_mem_hdl);

  reg_buf->transport_id = trans_hdl->id;
  reg_buf->buffer_owner = trans_hdl->me;
  reg_buf->ops = ops;
  reg_buf->payload_size = size;
  reg_buf->payload = (uint64_t) buffer;
  reg_buf->transport_private = (uint64_t) bgpdcmf_mem_hdl;

  if (peer != NULL)
    {
      reg_buf->buffer_owner = *peer;
    }

/*
  log_debug (nnti_debug_level, "rpc_buffer->payload_size=%ld",
	     reg_buf->payload_size);
*/
  reg_buf->buffer_addr.transport_id = NNTI_TRANSPORT_DCMF;
  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.size = size;
  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.buf = (uint64_t) buffer;
  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.owner_rank =
    DCMF_Messager_rank ();
  if (bgpdcmf_mem_hdl->ref_count == 1)
    {
      if (ops == NNTI_RECV_QUEUE)
	{
	  bgpdcmf_request_queue_handle *
	    q_hdl = &transport_global_data.req_queue;

	  /*
	   * This is a receive-only buffer.  This buffer is divisible by
	   * NNTI_REQUEST_BUFFER_SIZE.  This buffer can hold more than
	   * one short request.  Assume this buffer is a request queue.
	   */
	  bgpdcmf_mem_hdl->type = REQUEST_BUFFER;
	  memset (q_hdl, 0, sizeof (bgpdcmf_request_queue_handle));
	  q_hdl->reg_buf = reg_buf;	/* req buffer can be accessed from global req_queue  */
	  transport_global_data.req_queue.reg_buf = reg_buf;
	  //log_debug(nnti_debug_level, "register req buffer for server");

	  server_req_queue_init (q_hdl, buffer, size, num_elements);

	  reg_buf->payload_size = q_hdl->req_size;
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.type =
	    NNTI_DCMF_REQUEST_BUFFER;

	}
      else if (ops == NNTI_RECV_DST)
	{
	  if (size == NNTI_RESULT_BUFFER_SIZE)
	    {
	      /*
	       * This is a receive-only buffer.  This buffer can hold exactly
	       * one short result.  Assume this buffer is a result queue.
	       */
	      bgpdcmf_mem_hdl->type = RESULT_BUFFER;

	      rc =
		register_memory (bgpdcmf_mem_hdl, buffer,
				 NNTI_RESULT_BUFFER_SIZE, ops);

	      reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.type =
		NNTI_DCMF_RECEIVE_DST;
	    }
	  else
	    {
	      /*
	       * This is a receive-only buffer.  This buffer doesn't look
	       * like a request buffer or a result buffer.  I don't know
	       * what it is.  Assume it is a regular data buffer.
	       */
	      bgpdcmf_mem_hdl->type = RECEIVE_BUFFER;

	      rc = register_memory (bgpdcmf_mem_hdl, buffer, size, ops);
	      reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.type =
		NNTI_DCMF_RECEIVE_DST;

	    }

	}
      else if (ops == NNTI_SEND_SRC)
	{
	  bgpdcmf_mem_hdl->type = SEND_BUFFER;

	  rc = register_memory (bgpdcmf_mem_hdl, buffer, size, ops);
	  if (rc != NNTI_OK)
	    {
	      fprintf (stderr,
		       "failed registering short request in register memory\n");
	    }

	}
      else if (ops == NNTI_GET_DST)
	{
	  bgpdcmf_mem_hdl->type = GET_DST_BUFFER;
	  rc = register_memory (bgpdcmf_mem_hdl, buffer, size, ops);

	}
      else if (ops == NNTI_GET_SRC)
	{
	  bgpdcmf_mem_hdl->type = GET_SRC_BUFFER;


	  rc = register_memory (bgpdcmf_mem_hdl, buffer, size, ops);

	}
      else if (ops == NNTI_PUT_SRC)
	{
	  bgpdcmf_mem_hdl->type = PUT_SRC_BUFFER;


	  rc = register_memory (bgpdcmf_mem_hdl, buffer, size, ops);
	}
      else if (ops == NNTI_PUT_DST)
	{
	  bgpdcmf_mem_hdl->type = PUT_DST_BUFFER;

	  rc = register_memory (bgpdcmf_mem_hdl, buffer, size, ops);

	}
      else if (ops == (NNTI_GET_SRC | NNTI_PUT_DST))
	{
	  bgpdcmf_mem_hdl->type = RDMA_TARGET_BUFFER;

	  rc = register_memory (bgpdcmf_mem_hdl, buffer, size, ops);
	}
      else
	{
	  bgpdcmf_mem_hdl->type = UNKNOWN_BUFFER;
	}
    }
  wr = (dcmf_work_request *) calloc (1, sizeof (dcmf_work_request));
  wr->reg_buf = reg_buf;
  wr->last_op = DCMF_OP_REGISTER_RDMA;
  register_wc (wr);

  nthread_lock (&bgpdcmf_mem_hdl->wr_queue_lock);
  bgpdcmf_mem_hdl->wr_queue.push_back (wr);
  nthread_unlock (&bgpdcmf_mem_hdl->wr_queue_lock);

  if (rc == NNTI_OK)
    {
      reg_buf->buffer_addr.transport_id = NNTI_TRANSPORT_DCMF;
      memcpy (&reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.mem_hdl,
	      &bgpdcmf_mem_hdl->mem_hdl, sizeof (DCMF_Memregion_t));
      if (bgpdcmf_mem_hdl->type == REQUEST_BUFFER)
	{
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.size =
	    transport_global_data.req_queue.req_size;
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.buf =
	    (uint64_t) transport_global_data.req_queue.req_buffer;
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.type =
	    NNTI_DCMF_REQUEST_BUFFER;
	}
      else
	{
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.size =
	    reg_buf->payload_size;
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.buf =
	    (uint64_t) reg_buf->payload;
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.type =
	    NNTI_DCMF_SEND_SRC;

	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.wc_addr =
	    (uint64_t) & wr->wc;
	  memcpy (&reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.
		  wc_mem_hdl, &wr->wc_mem_hdl, sizeof (DCMF_Memregion_t));
	}

    }

  if (bgpdcmf_mem_hdl->ref_count == 1)
    {
      insert_buf_bufhash (reg_buf);
      //log_debug (nnti_debug_level, "bgpdcmf_mem_hdl->type==%llu",
		 //(uint64_t) bgpdcmf_mem_hdl->type);
      //log_debug (nnti_debug_level, "reg_buf.buf.hash==%llu",
		 //(uint64_t) hash6432shift (reg_buf->buffer_addr.
		//			   NNTI_remote_addr_t_u.bgpdcmf.buf));
    }

  return (rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t
NNTI_bgpdcmf_unregister_memory (NNTI_buffer_t * reg_buf)
{
  NNTI_result_t
    rc = NNTI_OK;
  bgpdcmf_memory_handle *
    bgpdcmf_mem_hdl = NULL;

  assert (reg_buf);

  //log_debug (nnti_debug_level, "enter unregistered buffer %p", reg_buf);

  bgpdcmf_mem_hdl = (bgpdcmf_memory_handle *) reg_buf->transport_private;

  assert (bgpdcmf_mem_hdl);
  bgpdcmf_mem_hdl->ref_count--;

  if (bgpdcmf_mem_hdl->ref_count == 0)
    {
      //log_debug (nnti_debug_level,
		 //"bgpdcmf_mem_hdl->ref_count is 0.  release all resources.");

      if (bgpdcmf_mem_hdl->type == REQUEST_BUFFER)
	{
	  server_req_queue_destroy (&transport_global_data.req_queue);

	}
      else
	{
	  unregister_memory (bgpdcmf_mem_hdl);
	}

      del_buf_bufhash (reg_buf);
      nthread_lock (&bgpdcmf_mem_hdl->wr_queue_lock);
      while (!bgpdcmf_mem_hdl->wr_queue.empty ())
	{
	  dcmf_work_request *
	    wr = bgpdcmf_mem_hdl->wr_queue.front ();
	  //log_debug (nnti_debug_level, "removing pending wr=%p", wr);
	  bgpdcmf_mem_hdl->wr_queue.pop_front ();
	  del_wr_wrhash (wr);
  	  //log_debug (nnti_debug_level, " called in unregister memory");
	}
      nthread_unlock (&bgpdcmf_mem_hdl->wr_queue_lock);
    }
  nthread_lock_fini (&bgpdcmf_mem_hdl->wr_queue_lock);
  if (bgpdcmf_mem_hdl)
    delete
      bgpdcmf_mem_hdl;

  reg_buf->transport_id = NNTI_TRANSPORT_NULL;
  reg_buf->payload_size = 0;
  reg_buf->payload = 0;
  reg_buf->transport_private = 0;

  //log_debug (nnti_debug_level, "exit");

  return (rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * Put the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t
NNTI_bgpdcmf_put (const NNTI_buffer_t * src_buffer_hdl,
		  const uint64_t src_offset,
		  const uint64_t src_length,
		  const NNTI_buffer_t * dest_buffer_hdl,
		  const uint64_t dest_offset)
{
  bgpdcmf_memory_handle *
    src_buf_mem_hdl = NULL;
  DCMF_Protocol_t
    put_prot;
  NNTI_result_t
    rc = NNTI_OK;
  dcmf_work_request *
    wr = NULL;

  wr = (dcmf_work_request *) calloc (1, sizeof (dcmf_work_request));
  assert (wr);

  wr->reg_buf = src_buffer_hdl;

  src_buf_mem_hdl =
    (bgpdcmf_memory_handle *) src_buffer_hdl->transport_private;
  size_t
    srv_rank =
    dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.owner_rank;
  {
    DCMF_Put_Configuration_t
    put_conf = { DCMF_DEFAULT_PUT_PROTOCOL };
    DCMF_CriticalSection_enter (0);
    DCMF_Put_register (&put_prot, &put_conf);
    DCMF_CriticalSection_exit (0);
  }
  wr->peer_rank = srv_rank;
  register_wc (wr);

  wr->op_state = RDMA_WRITE_INIT;
  wr->wc.op = DCMF_OP_PUT_TARGET;
  wr->wc.req_rank = transport_global_data.myrank;
  wr->wc.byte_len = src_length;
  wr->wc.src_offset = src_offset;
  wr->wc.dest_offset = dest_offset;
  wr->last_op = DCMF_OP_PUT_INITIATOR;

  wr->wc.local_send_complete = 1;
  wr->wc.remote_put_complete = 1;
  DCMF_Callback_t
  cb0 = { decrement, (void *) &wr->wc.local_send_complete };
  DCMF_Callback_t
    _cb_done;
  DCMF_Request_t
    req0;

  _cb_done.function = decrement;
  _cb_done.clientdata = (void *) &wr->wc.remote_put_complete;

/*
  size_t  length;
  void * base;
   DCMF_Memregion_query (&src_buf_mem_hdl->mem_hdl, &length, &base);
   fprintf(stderr, "bgpdcmf put SRC base %p size %ld\n", base, length);
  DCMF_Memregion_query ((DCMF_Memregion_t *) &dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.mem_hdl, &length, &base);
   fprintf(stderr, "bgpdcmf put DEST base %p size %ld\n", base, length);
*/

  nthread_lock (&nnti_bgpdcmf_lock);

  DCMF_CriticalSection_enter (0);
  DCMF_Put (&put_prot, &req0, cb0, DCMF_SEQUENTIAL_CONSISTENCY,
	    srv_rank, src_length, &src_buf_mem_hdl->mem_hdl,
	    (DCMF_Memregion_t *) & dest_buffer_hdl->buffer_addr.
	    NNTI_remote_addr_t_u.bgpdcmf.mem_hdl, src_offset, dest_offset,
	    _cb_done);
  while (wr->wc.remote_put_complete)
    DCMF_Messager_advance ();
  DCMF_CriticalSection_exit (0);
  //log_debug(nnti_debug_level, "messeger advance in put");
  memcpy (&wr->wc_dest_mem_hdl,
	  &dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.
	  wc_mem_hdl, sizeof (DCMF_Memregion_t));
  nthread_unlock (&nnti_bgpdcmf_lock);

  nthread_lock (&src_buf_mem_hdl->wr_queue_lock);
  src_buf_mem_hdl->wr_queue.push_back (wr);
  nthread_unlock (&src_buf_mem_hdl->wr_queue_lock);
  insert_wr_wrhash (wr);
  //log_debug (nnti_debug_level, "exit");

  return (rc);
}


/**
 * @brief Send a message to a peer.
 *
 * Send a message (<tt>msg_hdl</tt>) to a peer (<tt>peer_hdl</tt>).  It is expected that the
 * message is small, but the exact maximum size is transport dependent.
 */
NNTI_result_t
NNTI_bgpdcmf_send (const NNTI_peer_t * peer_hdl,
		   const NNTI_buffer_t * msg_hdl,
		   const NNTI_buffer_t * dest_hdl)
{
  NNTI_result_t
    rc = NNTI_OK;


  bgpdcmf_memory_handle *
    bgpdcmf_mem_hdl = NULL;


  assert (peer_hdl);
  assert (msg_hdl);

  bgpdcmf_mem_hdl = (bgpdcmf_memory_handle *) msg_hdl->transport_private;

  if ((dest_hdl == NULL) || (dest_hdl->ops == NNTI_RECV_QUEUE))
    {
      bgpdcmf_connection *
	conn =
	get_conn_rank (peer_hdl->peer.NNTI_remote_process_t_u.bgpdcmf.
		       pset_rank);
      assert (conn);
      request_send (&conn->queue_local_attrs,
		    &conn->queue_remote_attrs.server, msg_hdl,
		    conn->peer_rank);
      log_debug (nnti_debug_level, "sending request to (%s)", peer_hdl->url);

    }
  else
    {
      rc = NNTI_bgpdcmf_put (msg_hdl, 0, msg_hdl->payload_size, dest_hdl, 0);
      if (rc != NNTI_OK)
	log_error (nnti_debug_level, "Put() failed: %d", rc);
    }


  //log_debug (nnti_debug_level, "exit");

  return (rc);
}


/**
 * @brief Transfer data from a peer.
 *
 * Get the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t
NNTI_bgpdcmf_get (const NNTI_buffer_t * src_buffer_hdl,
		  const uint64_t src_offset,
		  const uint64_t src_length,
		  const NNTI_buffer_t * dest_buffer_hdl,
		  const uint64_t dest_offset)
{
  DCMF_Protocol_t
    get_protocol;
  NNTI_result_t
    rc = NNTI_OK;
  bgpdcmf_memory_handle *
    target_buf_mem_hdl;
  int
    srv_rank;
  dcmf_work_request *
    wr = NULL;


  wr = (dcmf_work_request *) calloc (1, sizeof (dcmf_work_request));
  wr->reg_buf = dest_buffer_hdl;
  srv_rank =
    src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.owner_rank;

  target_buf_mem_hdl =
    (bgpdcmf_memory_handle *) dest_buffer_hdl->transport_private;
  {
    DCMF_Result
      dcmf_result;
    DCMF_Get_Configuration_t
      conf;

    conf.protocol = DCMF_DEFAULT_GET_PROTOCOL;
    conf.network = DCMF_TORUS_NETWORK;

    dcmf_result = DCMF_Get_register (&get_protocol, &conf);
    assert (dcmf_result == DCMF_SUCCESS);
  }
/*
  DCMF_Memregion_t my_memregion;
  DCMF_CriticalSection_enter(0);
  DCMF_Result dcmf_result = DCMF_Memregion_create( &my_memregion, &bytes_out, src_length, dest_buffer_hdl->payload, 0 );
  DCMF_CriticalSection_exit(0);
  assert( dcmf_result==DCMF_SUCCESS && bytes_out==src_length );
*/
  wr->op_state = RDMA_READ_INIT;
  wr->wc.op = DCMF_OP_GET_TARGET;
  wr->wc.req_rank = srv_rank;
  wr->wc.byte_len = src_length;
  wr->wc.src_offset = src_offset;
  wr->wc.dest_offset = dest_offset;
  wr->last_op = DCMF_OP_GET_INITIATOR;
  wr->wc.local_get_complete = 1;
  wr->peer_rank = srv_rank;
  register_wc (wr);
  DCMF_Request_t
    request;
  DCMF_Callback_t
    done_callback;
  DCMF_Result
    dcmf_result;

  size_t
      regsize;
    void *
      base;
    DCMF_Memregion_query ((DCMF_Memregion_t *)&src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.mem_hdl, &regsize,
                          &base);
    log_debug(nnti_debug_level, "GET SRC base %p size %ld rank %d\n", base, regsize, srv_rank);
    DCMF_Memregion_query (&target_buf_mem_hdl->mem_hdl, &regsize, &base);
    log_debug(nnti_debug_level, "GET DEST base %p size %ld length %llu\n", base, regsize, src_length);

  nthread_lock (&nnti_bgpdcmf_lock);
  DCMF_CriticalSection_enter (0);
  done_callback.function = decrement;
  done_callback.clientdata = (void *) &wr->wc.local_get_complete;
  dcmf_result = DCMF_Get (&get_protocol,
			  &request,
			  done_callback,
			  DCMF_SEQUENTIAL_CONSISTENCY,
			  srv_rank,
			  src_length,
			  (DCMF_Memregion_t *) & src_buffer_hdl->buffer_addr.
			  NNTI_remote_addr_t_u.bgpdcmf.mem_hdl,
			  &target_buf_mem_hdl->mem_hdl, src_offset,
			  dest_offset);
  while (wr->wc.local_get_complete > 0)
    DCMF_Messager_advance ();
  DCMF_CriticalSection_exit (0);
  log_debug(nnti_debug_level, "messager advance in get");
  assert (dcmf_result == DCMF_SUCCESS);
  memcpy (&wr->wc_dest_mem_hdl,
	  &src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.
	  wc_mem_hdl, sizeof (DCMF_Memregion_t));
  nthread_unlock (&nnti_bgpdcmf_lock);
  nthread_lock (&target_buf_mem_hdl->wr_queue_lock);
  target_buf_mem_hdl->wr_queue.push_back (wr);
  nthread_unlock (&target_buf_mem_hdl->wr_queue_lock);
  insert_wr_wrhash (wr);
  //log_debug (nnti_debug_level, "exit");

  return (rc);

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
NNTI_result_t
NNTI_bgpdcmf_wait (const NNTI_buffer_t * reg_buf,
		   const NNTI_buf_ops_t remote_op,
		   const int timeout, NNTI_status_t * status)
{
  NNTI_result_t
    nnti_rc = NNTI_OK;
  bgpdcmf_memory_handle *
    bgpdcmf_mem_hdl = NULL;
  dcmf_work_request *
    wr = NULL;
  bgpdcmf_connection *
    conn = NULL;
  int
    timeout_per_call;
  uint8_t
    retry_count = 0;



  bgpdcmf_mem_hdl = (bgpdcmf_memory_handle *) reg_buf->transport_private;

  if (timeout < 0)
    timeout_per_call = MIN_TIMEOUT;
  else
    timeout_per_call = (timeout < MIN_TIMEOUT) ? MIN_TIMEOUT : timeout;

  retry_count = 0;
  wr = bgpdcmf_mem_hdl->wr_queue.front ();
  bgpdcmf_mem_hdl->wr_queue.pop_front ();
  if ((wr->last_op == DCMF_OP_REGISTER_RDMA))
    {
      if ((remote_op == NNTI_SEND_SRC) || (remote_op == NNTI_GET_DST)
	  || (remote_op == NNTI_PUT_SRC))
	{
	  wr = bgpdcmf_mem_hdl->wr_queue.front ();
	  bgpdcmf_mem_hdl->wr_queue.pop_front ();
	}
    }
  log_debug (nnti_debug_level, "bgpdcmf_wait  buffer = %p work request =%p", reg_buf, wr);
  nnti_rc = process_event (reg_buf, remote_op, wr, timeout_per_call);
  conn = get_conn_rank (wr->peer_rank);
  status->op = remote_op;
  status->result = nnti_rc;
  if (nnti_rc == NNTI_OK)
    {
      status->start = (uint64_t) reg_buf->payload;
      status->offset = wr->wc.byte_offset;
      status->length = wr->wc.byte_len;
      switch (remote_op)
	{
	case NNTI_SEND_SRC:
	case NNTI_GET_SRC:
	case NNTI_PUT_SRC:	/* I am client here */
	  create_peer (&status->src, transport_global_data.myrank);
	  create_peer (&status->dest, conn->peer_rank);
	  break;
	case NNTI_RECV_QUEUE:
	case NNTI_RECV_DST:
	case NNTI_GET_DST:
	case NNTI_PUT_DST:
	  create_peer (&status->src, conn->peer_rank);
	  create_peer (&status->dest, transport_global_data.myrank);
	  break;

	}
    }
  if (nnti_rc == NNTI_OK)
    {

      if ((bgpdcmf_mem_hdl->type == RDMA_TARGET_BUFFER)
	  || (bgpdcmf_mem_hdl->type == RECEIVE_BUFFER)
	  || (bgpdcmf_mem_hdl->type == RESULT_BUFFER)
	  || (bgpdcmf_mem_hdl->type == GET_SRC_BUFFER)
	  || (bgpdcmf_mem_hdl->type == REQUEST_BUFFER)
	  || (bgpdcmf_mem_hdl->type == PUT_DST_BUFFER))
	{
	  repost_recv_work_request((NNTI_buffer_t *)reg_buf, wr);
	}
      else
	{
	  del_wr_wrhash (wr);
	  log_debug(nnti_debug_level, " delete wr called in bgpdcmf_wait wr = %p buffer = %p", wr, reg_buf);
	  free (wr);
	}
    }


  return (nnti_rc);
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
NNTI_result_t
NNTI_bgpdcmf_waitany (const NNTI_buffer_t ** buf_list,
		      const uint32_t buf_count,
		      const NNTI_buf_ops_t remote_op,
		      const int timeout,
		      uint32_t * which, NNTI_status_t * status)
{
  NNTI_result_t
    nnti_rc = NNTI_OK;

  log_level
    debug_level = nnti_debug_level;

  log_debug (debug_level, "enter");

  assert (buf_list);
  assert (buf_count > 0);
  if (buf_count > 1)
    {
      /* if there is more than 1 buffer in the list, none of them can be a REQUEST_BUFFER */
      for (uint32_t i = 0; i < buf_count; i++)
	{
	  if (buf_list[i] != NULL)
	    {
	      assert (((bgpdcmf_memory_handle *) buf_list[i]->
		       transport_private)->type != REQUEST_BUFFER);
	    }
	}
    }
  assert (status);

  if (buf_count == 1)
    {
      nnti_rc = NNTI_bgpdcmf_wait (buf_list[0], remote_op, timeout, status);
      *which = 0;
      goto cleanup;
    }

  while (1)
    {
      if (is_any_buf_op_complete (buf_list, buf_count, which) == TRUE)
	{
	  log_debug (debug_level,
		     "buffer op already complete (which=%u, buf_list[%d]=%p)",
		     *which, *which, buf_list[*which]);
	  nnti_rc = NNTI_OK;
	}
      else
	{
	  log_debug (debug_level, "buffer op NOT complete (buf_list=%p)",
		     buf_list);

	  for (uint32_t i = 0; i < buf_count; i++)
	    {
	      if (buf_list[i] != NULL)
		{
		  nnti_rc =
		    NNTI_bgpdcmf_wait (buf_list[i], remote_op,
				       timeout / buf_count, status);
		  if (nnti_rc == NNTI_OK)
		    {
		      *which = i;
		      break;
		    }
		}
	    }

	}
    }

cleanup:

  //log_debug (debug_level, "exit");

  trios_stop_timer ("NNTI_bgpdcmf_waitany", total_time);

  return (nnti_rc);
}

NNTI_result_t
NNTI_bgpdcmf_waitall (const NNTI_buffer_t ** buf_list,
		      const uint32_t buf_count,
		      const NNTI_buf_ops_t remote_op,
		      const int timeout, NNTI_status_t ** status)
{
  return NNTI_OK;
}

/**
 * @brief Disable this transport.
 *
 * Shutdown the transport.  Any outstanding sends, gets and puts will be
 * canceled.  Any new transport requests will fail.
 *
 */
/**
 * @brief Disable this transport.
 *
 * Shutdown the transport.  Any outstanding sends, gets and puts will be
 * canceled.  Any new transport requests will fail.
 *
 */
NNTI_result_t
NNTI_bgpdcmf_fini (const NNTI_transport_t * trans_hdl)
{
  /* Free up injection counter and global inj fifo and reception FIFO */
  DCMF_Messager_finalize ();
  return (NNTI_OK);
}

static NNTI_result_t
register_wc (dcmf_work_request * wr)
{
  DCMF_Result
    dcmf_result;
  size_t
    bytes_out;


/* This is for Work completion memory handle where ACK is received after PUT/GET */
  dcmf_result =
    DCMF_Memregion_create ((DCMF_Memregion_t *) & wr->wc_mem_hdl, &bytes_out,
			   sizeof (nnti_bgpdcmf_work_completion), &wr->wc, 0);

  if (dcmf_result != DCMF_SUCCESS)
    {
      fprintf (stderr,
	       "DCMF memregion create failed in register memory Work completion handle\n");
      return (NNTI_EIO);
    }

  //log_debug (nnti_debug_level, "exit  wr(%p)", wr);

  return ((NNTI_result_t) DCMF_SUCCESS);

}


static
  NNTI_result_t
register_memory (bgpdcmf_memory_handle * hdl, void *buf, uint64_t len,
		 const NNTI_buf_ops_t remote_op)
{
  NNTI_result_t
    rc = NNTI_OK;		/* return code */
  DCMF_Result
    dcmf_result;
  size_t
    bytes_out;

  assert (hdl);
  dcmf_result =
    DCMF_Memregion_create (&hdl->mem_hdl, &bytes_out, len, buf, 0);
  if (dcmf_result != DCMF_SUCCESS)
    fprintf (stderr, "DCMF memregion create failed in register memory \n");

  return (rc);
}

static int
unregister_memory (bgpdcmf_memory_handle * hdl)
{
  int
    rc = NNTI_OK;		/* return code */

  log_debug (nnti_debug_level, "exit  hdl(%p)", hdl);

  return (rc);
}

static
  NNTI_result_t
send_ack (const NNTI_buffer_t * reg_buf, dcmf_work_request * wr)
{
  NNTI_result_t
    rc = NNTI_OK;
  bgpdcmf_memory_handle *
    bgpdcmf_mem_hdl = NULL;
  DCMF_Protocol_t
    put_prot;


  bgpdcmf_mem_hdl = (bgpdcmf_memory_handle *) reg_buf->transport_private;
#ifndef HAS_MPI
  DCMF_Messager_initialize ();
#endif
  {				/* init  put  */
    DCMF_Put_Configuration_t
    put_conf = { DCMF_DEFAULT_PUT_PROTOCOL };
    DCMF_CriticalSection_enter (0);
    DCMF_Put_register (&put_prot, &put_conf);
    DCMF_CriticalSection_exit (0);
  }
  nthread_lock(&nnti_bgpdcmf_lock);
  DCMF_CriticalSection_enter (0);
  volatile unsigned
    active0 = 1;
  DCMF_Callback_t
  cb0 = { decrement, (void *) &active0 };
  DCMF_Callback_t
    _cb_null;
  DCMF_Request_t
    req0;

  _cb_null.function = NULL;
  _cb_null.clientdata = (void *) NULL;

  /*  SHYAMALI 
     size_t  length;
     void * base;
     DCMF_Memregion_query (&wr->wc_mem_hdl, &length, &base);
     DCMF_Memregion_query ((DCMF_Memregion_t *) &wr->wc_dest_mem_hdl, &length, &base);
   */

  DCMF_Put (&put_prot, &req0, cb0, DCMF_SEQUENTIAL_CONSISTENCY,
	    wr->peer_rank, sizeof (nnti_bgpdcmf_work_completion),
	    &wr->wc_mem_hdl, &wr->wc_dest_mem_hdl, 0, 0, _cb_null);
  while (active0)
    DCMF_Messager_advance ();
  DCMF_CriticalSection_exit (0);
  nthread_unlock(&nnti_bgpdcmf_lock);
  //log_debug(nnti_debug_level, "messager advance in send ack");
  return (rc);
}


static
  NNTI_result_t
process_event (const NNTI_buffer_t * reg_buf,
	       const NNTI_buf_ops_t remote_op,
	       dcmf_work_request * wr, const int timeout)
{
  NNTI_result_t
    nnti_rc = NNTI_OK;
  bgpdcmf_memory_handle *
    bgpdcmf_mem_hdl = NULL;

  log_level
    debug_level = nnti_debug_level;

  bgpdcmf_mem_hdl = (bgpdcmf_memory_handle *) reg_buf->transport_private;

  //log_debug (nnti_debug_level, "enter");
  long
    entry_time = trios_get_time_ms ();
  long
    elapsed_time = 0;


  debug_level = nnti_debug_level;
  switch (bgpdcmf_mem_hdl->type)
    {
    case SEND_BUFFER:
      while (wr->wc.remote_put_complete > 0)
	{
          usleep(100);	
          log_debug(nnti_debug_level, "messager advance in process event"); 
	  elapsed_time = (trios_get_time_ms () - entry_time);
	  if (((timeout >= 0) && (elapsed_time >= timeout))
	      || trios_exit_now ())
	    {
	      log_debug (nnti_debug_level,
			 "timed out SEND BUFFER...timeout(%d) elapsed_time(%d) exit_now(%d)",
			 timeout, elapsed_time, trios_exit_now ());
	      nnti_rc = NNTI_ETIMEDOUT;
	      break;
	    }
	}
      wr->op_state = SEND_COMPLETE;
      wr->wc.op_complete = 1;
      if (wr->wc.op != DCMF_OP_SEND_REQ)
	{
	  fprintf (stderr, "process_event send buffer sending ack\n");
	  send_ack (reg_buf, wr);
	}
       log_debug (nnti_debug_level, "process event send buffer");
      break;
    case PUT_SRC_BUFFER:
      wr->last_op = DCMF_OP_PUT_INITIATOR;
      if (wr->op_state == RDMA_WRITE_INIT)
	{
	  while (wr->wc.remote_put_complete > 0)
	    {
	      elapsed_time = (trios_get_time_ms () - entry_time);
	      if (((timeout >= 0) && (elapsed_time >= timeout))
		  || trios_exit_now ())
		{
		  log_debug (nnti_debug_level,
			     "timed out PUT SRC BUFFER...timeout(%d) elapsed_time(%d) exit_now(%d)",
			     timeout, elapsed_time, trios_exit_now ());
		  nnti_rc = NNTI_ETIMEDOUT;
		  break;
		}

	    }
	  wr->wc.op_complete = 1;
	  wr->op_state = RDMA_COMPLETE;
          log_debug (nnti_debug_level, "process event put src  buffer");
	  send_ack (reg_buf, wr);
	}
      break;
    case GET_DST_BUFFER:
      wr->last_op = DCMF_OP_GET_INITIATOR;
      if (wr->op_state == RDMA_READ_INIT)
	{
	  while (wr->wc.local_get_complete > 0)
	    {
	      elapsed_time = (trios_get_time_ms () - entry_time);
	      if (((timeout >= 0) && (elapsed_time >= timeout))
		  || trios_exit_now ())
		{
		  log_debug (nnti_debug_level,
			     "timed out GET DST BUFFER...timeout(%d) elapsed_time(%d) exit_now(%d)",
			     timeout, elapsed_time, trios_exit_now ());
		  nnti_rc = NNTI_ETIMEDOUT;
		  break;
		}

	    }
	  wr->op_state = RDMA_COMPLETE;
	  wr->wc.op_complete = 1;
          log_debug (nnti_debug_level, "process event get dest buffer");
	  send_ack (reg_buf, wr);
	}
      break;
    case REQUEST_BUFFER:
      {
	uint64_t
	  index = 0;
	bgpdcmf_request_queue_handle *
	  q = &transport_global_data.req_queue;

	wr->last_op = DCMF_OP_NEW_REQUEST;

	log_debug(nnti_debug_level,
		   "recv completion - reg_buf=%p current_req_index =%llu processing=%llu\n",
		   reg_buf, q->req_index, q->req_processed);
	while (q->req_index == q->req_processed)
	  {
	 /*   elapsed_time = (trios_get_time_ms () - entry_time);
	    if (((timeout >= 0) && (elapsed_time >= timeout))
		|| trios_exit_now ())
	      {
		log_debug (nnti_debug_level,
			   "timed out REQUEST BUFFER...timeout(%d) elapsed_time(%d) exit_now(%d)",
			   timeout, elapsed_time, trios_exit_now ());
		nnti_rc = NNTI_ETIMEDOUT;
		break;
	      }
           */
		usleep(100);

	  }
	index = q->req_processed;
	/* Shyamali wait for work completion buffer to showup */
	while (q->wc_buffer[q->req_processed].req_received != 1)
	  {
	    usleep (100);
	  }

	q->wc_buffer[q->req_processed].byte_offset = index * q->req_size;
	q->wc_buffer[q->req_processed].byte_len = q->req_size;
	wr->wc = q->wc_buffer[q->req_processed];
	wr->peer_rank = q->wc_buffer[q->req_processed].req_rank;
	fprintf (stderr, "This request came from %d  process index %llu\n",
		 wr->peer_rank, q->req_processed);
	wr->op_state = RECV_COMPLETE;
	q->req_processed++;
	q->total_req_processed++;
	if (q->req_processed > q->req_count)
	  {
	    log_error (nnti_debug_level,
		       "req_processed(%llu) > req_count(%llu) fail",
		       q->req_processed, q->req_count);
	  }
	if (q->req_processed == (q->req_count / 2))
	  {
	    if (q->req_index >= q->req_count)
	      {
		reset_req_index (q);
		log_debug (nnti_debug_level,
			   "resetting req_processed(%llu) total_req_processed(%llu)",
			   q->req_processed, q->total_req_processed);
	      }
	    else
	      {
		log_debug (nnti_debug_level,
			   "skipping reset req_processed(%llu) total_req_processed(%llu)",
			   q->req_processed, q->total_req_processed);
	      }
	  }
	if (q->req_processed == q->req_processed_reset_limit)
	  {
	    q->req_processed = 0;
	  }
	log_debug (nnti_debug_level,
		   "current req_processed(%llu) req_count(%llu)",
		   q->req_processed, q->req_count);
      }
      break;
    case RECEIVE_BUFFER:
      wr->last_op = DCMF_OP_RECEIVE;
      log_debug (debug_level,
		 "receive buffer - recv completion - reg_buf==%p", reg_buf);

      while (wr->wc.op_complete != 1)
	{
	  elapsed_time = (trios_get_time_ms () - entry_time);
	  if (((timeout >= 0) && (elapsed_time >= timeout))
	      || trios_exit_now ())
	    {
	      log_debug (nnti_debug_level,
			 "timed out TARGET RECEIVE  BUFFER...timeout(%d) elapsed_time(%d) exit_now(%d)",
			 timeout, elapsed_time, trios_exit_now ());
	      nnti_rc = NNTI_ETIMEDOUT;
	      break;
	    }

	}
      wr->peer_rank = wr->wc.req_rank;
      wr->op_state = RECV_COMPLETE;
      break;
    case RESULT_BUFFER:
      wr->last_op = DCMF_OP_PUT_TARGET;
      while (wr->wc.op_complete != 1)
	{
	  elapsed_time = (trios_get_time_ms () - entry_time);
	  if (((timeout >= 0) && (elapsed_time >= timeout))
	      || trios_exit_now ())
	    {
	      log_debug (nnti_debug_level,
			 "timed out RESULT BUFFER...timeout(%d) elapsed_time(%d) exit_now(%d)",
			 timeout, elapsed_time, trios_exit_now ());
	      nnti_rc = NNTI_ETIMEDOUT;
	      break;
	    }

	}
      wr->peer_rank = wr->wc.req_rank;
      wr->op_state = RDMA_COMPLETE;
      log_debug(nnti_debug_level, "process event result buffer");
      break;
    case PUT_DST_BUFFER:
      wr->last_op = DCMF_OP_PUT_TARGET;
      //fprintf (stderr, "PUT DST Buffer\n");
      while (wr->wc.op_complete != 1)
	{
	  elapsed_time = (trios_get_time_ms () - entry_time);
	  if (((timeout >= 0) && (elapsed_time >= timeout))
	      || trios_exit_now ())
	    {
	      log_debug (nnti_debug_level,
			 "timed out PUT DST BUFFER...timeout(%d) elapsed_time(%d) exit_now(%d)",
			 timeout, elapsed_time, trios_exit_now ());
	      nnti_rc = NNTI_ETIMEDOUT;
	      break;
	    }

	}
      wr->op_state = RDMA_COMPLETE;
      log_debug(nnti_debug_level, "process event PUT DST  buffer");
      break;
    case GET_SRC_BUFFER:
      wr->last_op = DCMF_OP_GET_TARGET;
      //fprintf (stderr, "GET SRC BUFFER\n");
      if (wr->op_state == RDMA_READ_INIT)
	{
	  while (wr->wc.op_complete != 1)
	    {
	      elapsed_time = (trios_get_time_ms () - entry_time);
	      if (((timeout >= 0) && (elapsed_time >= timeout))
		  || trios_exit_now ())
		{
		  log_debug (nnti_debug_level,
			     "timed out GET SRC BUFFER...timeout(%d) elapsed_time(%d) exit_now(%d)",
			     timeout, elapsed_time, trios_exit_now ());
		  nnti_rc = NNTI_ETIMEDOUT;
		  break;
		}

	    }
	  wr->op_state = RDMA_COMPLETE;
	}
      log_debug(nnti_debug_level, "process event GET SRC  buffer");
      break;
    case RDMA_TARGET_BUFFER:
      if ((wr->last_op == DCMF_OP_GET_INITIATOR) ||
	  (wr->last_op == DCMF_OP_PUT_INITIATOR))
	{

	  if (wr->op_state == RDMA_TARGET_INIT)
	    {
	      while (wr->wc.op_complete != 1)
		{
		  elapsed_time = (trios_get_time_ms () - entry_time);
		  if (((timeout >= 0) && (elapsed_time >= timeout))
		      || trios_exit_now ())
		    {
		      log_debug (nnti_debug_level,
				 "timed out RDMA TARGET  BUFFER...timeout(%d) elapsed_time(%d) exit_now(%d)",
				 timeout, elapsed_time, trios_exit_now ());
		      nnti_rc = NNTI_ETIMEDOUT;
		      break;
		    }


		}

	      wr->op_state = RDMA_COMPLETE;

	    }
	}
      else
	{
	  if (wr->op_state == RDMA_TARGET_INIT)
	    {

	      log_debug (debug_level, "RDMA target completion - reg_buf==%p",
			 reg_buf);
	      while (wr->wc.op_complete != 1)
		{
		  elapsed_time = (trios_get_time_ms () - entry_time);
		  if (((timeout >= 0) && (elapsed_time >= timeout))
		      || trios_exit_now ())
		    {
		      log_debug (nnti_debug_level,
				 "timed out RDMA BUFFER...timeout(%d) elapsed_time(%d) exit_now(%d)",
				 timeout, elapsed_time, trios_exit_now ());
		      nnti_rc = NNTI_ETIMEDOUT;
		      break;
		    }

		}
	      wr->op_state = RDMA_COMPLETE;

	    }
	}

      break;
    case UNKNOWN_BUFFER:
      fprintf (stderr, "UNKNOWN_BUFFER\n");
      break;
    }
  //log_debug (nnti_debug_level, "exit");
  return (nnti_rc);
}


static void
create_peer (NNTI_peer_t * peer, int rank)
{
  sprintf (peer->url, "dcmf://%d", rank);

  peer->peer.transport_id = NNTI_TRANSPORT_DCMF;
  peer->peer.NNTI_remote_process_t_u.bgpdcmf.pset_rank = rank;
}


/*
 * Keep looping until all bytes have been accepted by the kernel.
 */
static
  NNTI_result_t
inject_full (int rank, const void *buf, size_t num)
{
  NNTI_result_t
    rc = NNTI_OK;
  DCQuad
    msginfo;
  DCMF_Request_t
    request;
  
  send_active = 1;
  nthread_lock(&nnti_bgpdcmf_lock);
  DCMF_CriticalSection_enter (0);
  DCMF_Callback_t
  cb_info = { decrement, (void *) &send_active };
  DCMF_Send (&send_prot,
	     &request,
	     cb_info,
	     DCMF_SEQUENTIAL_CONSISTENCY,
	     rank, num, (char *) buf, &msginfo, 1);

  TRACE_ERR ((stderr, "inject to  (%zd)  Before advance size = %d\n", rank,
	      num));
  while (send_active)
    DCMF_Messager_advance ();
  DCMF_CriticalSection_exit (0);
  nthread_unlock(&nnti_bgpdcmf_lock);
  send_active = 1;
  return rc;
}

static
  NNTI_result_t
new_client_connection (bgpdcmf_connection * c, int peer_rank)
{
  NNTI_result_t
    rc;

  /*
   * Values passed DMA MEMFIO for initial connection.
   */
  struct
  {
    nnti_bgpdcmf_server_queue_attrs
      server_attrs;
  } sa_in;
  struct
  {
    nnti_bgpdcmf_client_queue_attrs
      client_attrs;
  } ca_out;


  c->connection_type = CLIENT_CONNECTION;
  c->peer_rank = peer_rank;
  memset (&sa_in, 0, sizeof (sa_in));

  _recv_active = 1;
   confirm_conn = 1;
  nthread_lock(&nnti_bgpdcmf_lock);
  while (_recv_active)
    DCMF_Messager_advance ();
  nthread_unlock(&nnti_bgpdcmf_lock);
  memcpy (&sa_in, &_recv_buffer, sizeof (sa_in));
  _recv_active = 1;

  c->queue_remote_attrs.server = sa_in.server_attrs;
  client_req_queue_init (c);
  memset (&ca_out, 0, sizeof (ca_out));
  ca_out.client_attrs.req_index = c->queue_local_attrs.req_index;
  ca_out.client_attrs.req_index_addr = c->queue_local_attrs.req_index_addr;
  memcpy (&ca_out.client_attrs.req_index_mem_hdl,
	  &c->queue_local_attrs.req_index_mem_hdl, sizeof (DCMF_Memregion_t));

  rc = inject_full (peer_rank, &ca_out, sizeof (ca_out));
  nthread_lock(&nnti_bgpdcmf_lock);
  while (confirm_conn)
	DCMF_Messager_advance ();
  nthread_unlock(&nnti_bgpdcmf_lock);
  return rc;
}

static
  NNTI_result_t
new_server_connection (bgpdcmf_connection * c, int rank)
{
  NNTI_result_t
    rc;
  bgpdcmf_request_queue_handle *
    q_hdl = &transport_global_data.req_queue;

  /*
   * Values passed through DMA MEMFIFO interface to permit initial connection.
   */
  struct
  {
    nnti_bgpdcmf_server_queue_attrs
      server_attrs;
  } sa_out;

  assert (transport_global_data.req_queue.reg_buf);

  c->connection_type = SERVER_CONNECTION;
  c->peer_rank = rank;

  memset (&sa_out, 0, sizeof (sa_out));

  sa_out.server_attrs.req_index_addr =
    transport_global_data.req_queue.req_index_addr;

  memcpy (&sa_out.server_attrs.req_index_mem_hdl,
	  &transport_global_data.req_queue.req_index_mem_hdl,
	  sizeof (DCMF_Memregion_t));

  sa_out.server_attrs.req_buffer_addr =
    (uint64_t) transport_global_data.req_queue.req_buffer;
  sa_out.server_attrs.req_size = transport_global_data.req_queue.req_size;
  sa_out.server_attrs.req_count = transport_global_data.req_queue.req_count;
  memcpy (&sa_out.server_attrs.req_mem_hdl, q_hdl->req_mem_hdl,
	  sizeof (DCMF_Memregion_t));
  sa_out.server_attrs.wc_buffer_addr =
    (uint64_t) transport_global_data.req_queue.wc_buffer;
  memcpy (&sa_out.server_attrs.wc_mem_hdl,
	  &transport_global_data.req_queue.wc_mem_hdl,
	  sizeof (DCMF_Memregion_t));
  sa_out.server_attrs.unblock_buffer_addr =
    (uint64_t) transport_global_data.req_queue.unblock_buffer;
  memcpy (&sa_out.server_attrs.unblock_mem_hdl,
	  &transport_global_data.req_queue.unblock_mem_hdl,
	  sizeof (DCMF_Memregion_t));
  rc = inject_full (rank, &sa_out, sizeof (sa_out));
  if (rc)
    {
      fprintf (stderr, "new server connection: failed inject full\n");
      goto out;
    }

out:
  return rc;
}

static
  NNTI_result_t
insert_conn_rank (const uint32_t pset_rank, bgpdcmf_connection * conn)
{
  NNTI_result_t
    rc = NNTI_OK;

  nthread_lock (&nnti_conn_peer_lock);
  if (connections_by_rank.find (pset_rank) != connections_by_rank.end ())
    {
      if (connections_by_rank[pset_rank] == conn){
	  log_debug (nnti_debug_level, "connection already exists");
	  return(rc);
       }
    }
  connections_by_rank[pset_rank] = conn;	// add to connection map
  nthread_unlock (&nnti_conn_peer_lock);

  log_debug (nnti_debug_level, "peer connection added (conn=%p)", conn);

  return (rc);
}

static bgpdcmf_connection *
get_conn_rank (const uint32_t pset_rank)
{

  int
    key = pset_rank;
  bgpdcmf_connection *
    conn = NULL;

  nthread_lock (&nnti_conn_peer_lock);
  if (connections_by_rank.find (key) != connections_by_rank.end ())
    {
      conn = connections_by_rank[key];
    }

  nthread_unlock (&nnti_conn_peer_lock);

  if (conn != NULL)
    {
      return conn;
    }
  log_debug (nnti_debug_level, "connection NOT found");
  return (NULL);
}

static bgpdcmf_connection *
del_conn_rank (const uint32_t pset_rank)
{
  bgpdcmf_connection *
    conn = NULL;
  int
    key;


  key = pset_rank;

  nthread_lock (&nnti_conn_peer_lock);
  if (connections_by_rank.find (key) != connections_by_rank.end ())
    {
      conn = connections_by_rank[key];
    }

  nthread_unlock (&nnti_conn_peer_lock);

  if (conn != NULL)
    {
      connections_by_rank.erase (key);
    }
  else
    {
      log_debug (nnti_debug_level, "connection NOT found");
    }

  return (conn);
}

/**
 * @brief initialize
 */
static
  NNTI_result_t
init_connection (bgpdcmf_connection ** conn,
		 int peer_rank, const int is_server)
{
  NNTI_result_t
    rc = NNTI_OK;		/* return code */

  bgpdcmf_connection *
    c = NULL;


  c = (bgpdcmf_connection *) calloc (1, sizeof (bgpdcmf_connection));
  if (c == NULL)
    {
      log_error (nnti_debug_level,
		 "calloc returned NULL.  out of memory?: %s",
		 strerror (errno));
      rc = NNTI_ENOMEM;
      goto out;
    }

  if (is_server)
    {
      rc = new_server_connection (c, peer_rank);
    }
  else
    {
      rc = new_client_connection (c, peer_rank);
    }
  if (rc)
    {
      close_connection (c);
      c = NULL;
      goto out;
    }

  *conn = c;

out:
  return (rc);
}

/*
 * At an explicit BYE message, or at finalize time, shut down a connection.
 * If descriptors are posted, defer and clean up the connection structures
 * later.
 */
static void
close_connection (bgpdcmf_connection * c)
{

  if (c == NULL)
    return;

  log_debug (nnti_debug_level, "enter");

  if (c->connection_type == CLIENT_CONNECTION)
    {
      client_req_queue_destroy (c);
    }

  c->state = DISCONNECTED;

  log_debug (nnti_debug_level, "exit");
}

/**
 * Check for new connections.  The listening socket is left nonblocking
 * so this test can be quick; but accept is not really that quick compared
 * to polling an Gemini interface, for instance.  Returns >0 if an accept worked.
 */
static
  NNTI_result_t
check_poll_for_new_connections ()
{
  NNTI_result_t
    rc = NNTI_OK;

  bgpdcmf_connection *
    conn = NULL;
  uint32_t
    peer_rank;
  int
    i;

  _recv_active = 1;
   _is_server = 1;
  log_debug(nnti_debug_level, "server poll for new connection");
  while (_recv_active)
    {

  	nthread_lock (&nnti_bgpdcmf_lock);
  	for (i = 0; i< 1024; ++i)
    	{
        	peer_rank = transport_global_data.conn_req[i];
       		if(peer_rank != -1) {
       			conn = get_conn_rank(peer_rank);
       			if (conn == NULL) {
      				rc = init_connection (&conn, peer_rank, 1);
       				if (rc != NNTI_OK)
       				{
            				fprintf (stderr, "Server make new connection failed\n");
        			}
      				insert_conn_rank (peer_rank, conn);
       			}
		}
    	}
	  DCMF_Messager_advance ();
     	  nthread_unlock(&nnti_bgpdcmf_lock);
	  usleep(10000);
    }
cleanup:
  return rc;
}

/**
 * @brief Continually check for new connection attempts.
 *
 */
static void *
connection_listener_thread (void *args)
{
  NNTI_result_t
    rc = NNTI_OK;

  log_debug (nnti_debug_level,
	     "started thread to listen for client connection attempts");

  /* SIGINT (Ctrl-C) will get us out of this loop */
  while (1)
    {
      log_debug (nnti_debug_level, "listening for new connection");
      rc = check_poll_for_new_connections ();
      if (rc != NNTI_OK)
	{
	  log_fatal (nnti_debug_level,
		     "error returned from nssi_bgpdcmf_server_listen_for_client: %d",
		     rc);
	  continue;
	}
    }

  return (NULL);
}

static void *
fetch_req_index_thread (void *args)
{
  int
    i;
  size_t
    peer;
  DCMF_Protocol_t
    put_prot;
  bgpdcmf_request_queue_handle *
    q = &transport_global_data.req_queue;

  log_debug (nnti_debug_level, "started thread to service req_index fetch");
  {
    DCMF_Put_Configuration_t
    put_conf = { DCMF_DEFAULT_PUT_PROTOCOL };
    DCMF_CriticalSection_enter (0);
    DCMF_Put_register (&put_prot, &put_conf);
    DCMF_CriticalSection_exit (0);
  }
  volatile unsigned
    active0 = 1;
  DCMF_Callback_t
  cb0 = { decrement, (void *) &active0 };
  DCMF_Callback_t
    _cb_null;
  DCMF_Request_t
    req0;
  _cb_null.function = NULL;
  _cb_null.clientdata = (void *) NULL;


  /* SIGINT (Ctrl-C) will get us out of this loop */
  log_debug (nnti_debug_level, "listening for buffer index to send request");
  while (1)
    {
      for (i = 0; i < MAX_CONNECTION; ++i)
	{
	  if (q->unblock_buffer[i].req_rank == i)
	    {
	      nthread_lock (&nnti_index_lock);
	      peer = q->unblock_buffer[i].req_rank;
	      //fprintf (stderr, " fetch index rank %d  index %d\n", peer, i);
	      DCMF_CriticalSection_enter (0);
	      if (peer != -1)
		{
		  active0 = 1;

		  DCMF_Put (&put_prot, &req0, cb0,
			    DCMF_SEQUENTIAL_CONSISTENCY, peer,
			    sizeof (uint64_t), &q->req_index_mem_hdl,
			    (DCMF_Memregion_t *) & q->unblock_buffer[i].
			    response_hdl, 0, 0, _cb_null);
		  while (active0)
		    DCMF_Messager_advance ();
		  q->req_index++;
		  q->unblock_buffer[i].req_rank = -1;
		}
	      DCMF_CriticalSection_exit (0);
	      nthread_unlock (&nnti_index_lock);
	    }
	}
    }

  return (NULL);
}


/**
 * @brief Start a thread to check for new connection attempts.
 *
 */
static int
start_connection_listener_thread ()
{
  int
    rc = 0;
  pthread_t
    thread;

  /* Create the thread. Do we want special attributes for this? */
  rc = pthread_create (&thread, NULL, connection_listener_thread, NULL);
  if (rc)
    {
      log_error (nnti_debug_level,
		 "could not spawn thread for new connection");
      rc = NNTI_EBADRPC;
    }
  return rc;

}

static int
start_index_thread ()
{
  int
    rc = 0;
  pthread_t
    thread1;
  rc = pthread_create (&thread1, NULL, fetch_req_index_thread, NULL);
  if (rc)
    {
      log_error (nnti_debug_level,
		 "could not spawn thread for fetch_req buffer index");
      rc = NNTI_EBADRPC;
    }

  return rc;
}



static int
reset_req_index (bgpdcmf_request_queue_handle * req_queue_attrs)
{
  uint64_t
    value_before_reset = 0;

/* Do a remote get or Direct put to update req_index */

  req_queue_attrs->last_index_before_reset = value_before_reset;
  req_queue_attrs->req_processed_reset_limit = value_before_reset;
  req_queue_attrs->req_index = 0;
  log_debug (nnti_debug_level, "index before reset(%llu).",
	     req_queue_attrs->last_index_before_reset);
  log_debug (nnti_debug_level, "index after reset(%llu).",
	     req_queue_attrs->req_index);


  return (0);
}

static int
fetch_server_req_buffer_offset (nnti_bgpdcmf_client_queue_attrs *
				local_req_queue_attrs,
				nnti_bgpdcmf_server_queue_attrs *
				remote_req_queue_attrs, uint64_t addend,
				uint64_t * prev_offset, int srv_rank)
{


  uint32_t
    req_size = sizeof (nnti_bgpdcmf_fetch_index_req);
  DCMF_Protocol_t
    put_prot;
  {
    DCMF_Put_Configuration_t
    put_conf = { DCMF_DEFAULT_PUT_PROTOCOL };
    DCMF_CriticalSection_enter (0);
    DCMF_Put_register (&put_prot, &put_conf);
    DCMF_CriticalSection_exit (0);
  }

  local_req_queue_attrs->req_index = -1;
  volatile unsigned
    active0 = 1;
  DCMF_Callback_t
  cb0 = { decrement, (void *) &active0 };
  DCMF_Callback_t
    _cb_null;
  DCMF_Request_t
    req0;

  _cb_null.function = NULL;
  _cb_null.clientdata = (void *) NULL;

  size_t
    length;
  void *
    base;
  DCMF_Memregion_query (&local_req_queue_attrs->fetch_index_hdl, &length,
			&base);
  nthread_lock(&nnti_bgpdcmf_lock);
  DCMF_CriticalSection_enter (0);

  DCMF_Put (&put_prot, &req0, cb0, DCMF_SEQUENTIAL_CONSISTENCY,
	    srv_rank, sizeof (nnti_bgpdcmf_fetch_index_req),
	    &local_req_queue_attrs->fetch_index_hdl,
	    (DCMF_Memregion_t *) & remote_req_queue_attrs->unblock_mem_hdl, 0,
	    req_size * transport_global_data.myrank, _cb_null);
  while (active0)
    DCMF_Messager_advance ();
  DCMF_CriticalSection_exit (0);
  nthread_unlock(&nnti_bgpdcmf_lock);

  while (local_req_queue_attrs->req_index == -1)
    usleep (10);
  *prev_offset = local_req_queue_attrs->req_index;
  return (0);
}

static int
send_req (nnti_bgpdcmf_client_queue_attrs * local_req_queue_attrs,
	  nnti_bgpdcmf_server_queue_attrs * remote_req_queue_attrs,
	  uint64_t offset, const NNTI_buffer_t * reg_buf, int rank,
	  dcmf_work_request * wr)
{
  bgpdcmf_memory_handle *
    src_buf_mem_hdl = NULL;
  src_buf_mem_hdl = (bgpdcmf_memory_handle *) reg_buf->transport_private;
  uint64_t
    length;
  DCMF_Protocol_t
    put_prot;
  {
    DCMF_Put_Configuration_t
    put_conf = { DCMF_DEFAULT_PUT_PROTOCOL };
    DCMF_CriticalSection_enter (0);
    DCMF_Put_register (&put_prot, &put_conf);
    DCMF_CriticalSection_exit (0);
  }
  length = reg_buf->payload_size;
  {
    wr->wc.remote_put_complete = 1;	/* callback will decr  when recv done on remote side */
    wr->wc.local_send_complete = 1;	/* First callback returns after  injected into local TORUS */
    DCMF_Callback_t
    cb0 = { decrement, (void *) &wr->wc.local_send_complete };
    DCMF_Callback_t
      _cb_done;
    DCMF_Request_t
      req0;

    _cb_done.function = decrement;
    _cb_done.clientdata = (void *) &wr->wc.remote_put_complete;
/*
    size_t
      regsize;
    void *
      base;
    DCMF_Memregion_query (&remote_req_queue_attrs->req_mem_hdl, &regsize,
			  &base);
    log_debug(nnti_debug_level, "remote req send base %p size %ld\n", base, regsize);
    DCMF_Memregion_query (&src_buf_mem_hdl->mem_hdl, &regsize, &base);
    log_debug(nnti_debug_level, "local buffer req send base %p size %ld length %llu\n", base, regsize, length);
 
*/   
    nthread_lock(&nnti_bgpdcmf_lock);
    DCMF_CriticalSection_enter (0);
    DCMF_Put (&put_prot, &req0, cb0, DCMF_SEQUENTIAL_CONSISTENCY,
	      rank, length, &src_buf_mem_hdl->mem_hdl,
	      &remote_req_queue_attrs->req_mem_hdl, 0, offset, _cb_done);
    while (wr->wc.remote_put_complete)
      DCMF_Messager_advance ();
    DCMF_CriticalSection_exit (0);
    nthread_unlock(&nnti_bgpdcmf_lock);
  }
  return (0);
}

static int
send_req_wc (nnti_bgpdcmf_client_queue_attrs * client_q,
	     nnti_bgpdcmf_server_queue_attrs * server_q,
	     const NNTI_buffer_t * reg_buf, uint64_t offset,
	     dcmf_work_request * wr, int rank)
{
  bgpdcmf_memory_handle *
    bgpdcmf_mem_hdl = NULL;

  bgpdcmf_mem_hdl = (bgpdcmf_memory_handle *) reg_buf->transport_private;
  DCMF_Protocol_t
    put_prot;
  {
    DCMF_Put_Configuration_t
    put_conf = { DCMF_DEFAULT_PUT_PROTOCOL };
    DCMF_CriticalSection_enter (0);
    DCMF_Put_register (&put_prot, &put_conf);
    DCMF_CriticalSection_exit (0);
  }

  {
    volatile unsigned
      active0 = 1;
    DCMF_Callback_t
    cb0 = { decrement, (void *) &active0 };
    DCMF_Callback_t
      _cb_null;
    DCMF_Request_t
      req0;

    _cb_null.function = NULL;
    _cb_null.clientdata = (void *) NULL;
    register_wc (wr);
/*
  SHYAMALI DEBUG
  if(transport_global_data.myrank == 0) {
  	size_t  regsize;
  	void * base;
   	DCMF_Memregion_query (&server_q->wc_mem_hdl, &regsize, &base);
  	 fprintf(stderr, "server work completion base %p size %ld\n", base, regsize);
  	 DCMF_Memregion_query (&wr->wc_mem_hdl, &regsize, &base);
   	fprintf(stderr, "client  work completion base %p size %ld\n", base, regsize);
   }
 */
    nthread_lock(&nnti_bgpdcmf_lock);
    DCMF_CriticalSection_enter (0);
    DCMF_Put (&put_prot, &req0, cb0, DCMF_SEQUENTIAL_CONSISTENCY,
	      rank, sizeof (nnti_bgpdcmf_work_completion), &wr->wc_mem_hdl,
	      &server_q->wc_mem_hdl, 0, offset, _cb_null);
    while (active0)
      DCMF_Messager_advance ();
    DCMF_CriticalSection_exit (0);
    nthread_unlock(&nnti_bgpdcmf_lock);
  }
  return (0);

}



static int
request_send (nnti_bgpdcmf_client_queue_attrs * client_q,
	      nnti_bgpdcmf_server_queue_attrs * server_q,
	      const NNTI_buffer_t * reg_buf, int rank)
{
  int
    rc = 0;
  bgpdcmf_memory_handle *
    bgpdcmf_mem_hdl = NULL;
  uint64_t
    offset;
  dcmf_work_request *
    wr = NULL;

  uint32_t
    wc_size = sizeof (nnti_bgpdcmf_work_completion);

  bgpdcmf_mem_hdl = (bgpdcmf_memory_handle *) reg_buf->transport_private;

  assert (bgpdcmf_mem_hdl);
  wr = (dcmf_work_request *) malloc (sizeof (dcmf_work_request));
  memset (wr, 0, sizeof (dcmf_work_request));
  assert (wr);

  wr->reg_buf = reg_buf;

  bgpdcmf_mem_hdl->wr_queue.push_back (wr);

  insert_wr_wrhash (wr);
  rc = fetch_server_req_buffer_offset (client_q, server_q, 1, &offset, rank);

  rc =
    send_req (client_q, server_q, offset * server_q->req_size, reg_buf, rank,
	      wr);
  if (rc != NNTI_OK)
    log_error (nnti_debug_level, "send_req() failed: %d", rc);
  wr->wc.op = DCMF_OP_SEND_REQ;
  wr->wc.req_received = 1;
  wr->wc.req_rank = transport_global_data.myrank;
  wr->wc.byte_len = reg_buf->payload_size;
  wr->wc.byte_offset = client_q->req_index * server_q->req_size;
  wr->wc.src_offset = 0;
  wr->wc.dest_offset = client_q->req_index * server_q->req_size;
  wr->peer_rank = rank;

  log_debug (nnti_debug_level, "calling send_req_wc()");
  rc = send_req_wc (client_q, server_q, reg_buf, offset * wc_size, wr, rank);
  return (0);
}



static int
client_req_queue_init (bgpdcmf_connection * c)
{
  int
    rc = 0;
  size_t
    bytes;

  nnti_bgpdcmf_client_queue_attrs *
    q = &c->queue_local_attrs;
  q->req_index = 0;
  q->req_index_addr = (uint64_t) & q->req_index;
  DCMF_Memregion_create (&q->req_index_mem_hdl, &bytes,
			 sizeof (uint64_t), (void *) q->req_index_addr, 0);
  client_fetch_req.req_rank = transport_global_data.myrank;
  memcpy (&client_fetch_req.response_hdl, &q->req_index_mem_hdl,
	  sizeof (DCMF_Memregion_t));
  DCMF_Memregion_create (&q->fetch_index_hdl, &bytes,
			 sizeof (nnti_bgpdcmf_fetch_index_req),
			 (void *) &client_fetch_req, 0);
  return (rc);
}

static int
client_req_queue_destroy (bgpdcmf_connection * c)
{
  return (0);
}

static int
server_req_queue_init (bgpdcmf_request_queue_handle * q,
		       char *buffer, uint64_t req_size, uint64_t req_count)
{
  int
    i;
  size_t
    bytes;
  uint64_t
    unblock_size;
  bgpdcmf_memory_handle *
    bgpdcmf_mem_hdl = (bgpdcmf_memory_handle *) q->reg_buf->transport_private;


  q->req_buffer = buffer;
  q->req_size = req_size;
  q->req_count = req_count;
  q->req_buffer_size = req_count * req_size;

  q->wc_buffer_size = q->req_count * sizeof (nnti_bgpdcmf_work_completion);
  q->wc_buffer =
    (nnti_bgpdcmf_work_completion *) calloc (q->req_count,
					     sizeof
					     (nnti_bgpdcmf_work_completion));

  unblock_size = MAX_CONNECTION * sizeof (nnti_bgpdcmf_work_completion);
  q->unblock_buffer =
    (nnti_bgpdcmf_fetch_index_req *) calloc (MAX_CONNECTION,
					     sizeof
					     (nnti_bgpdcmf_fetch_index_req));
  DCMF_Memregion_create (&q->req_mem_hdl, &bytes, q->req_buffer_size, buffer,
			 0);
  DCMF_Memregion_create (&q->wc_mem_hdl, &bytes, q->wc_buffer_size,
			 q->wc_buffer, 0);
  DCMF_Memregion_create (&q->unblock_mem_hdl, &bytes, unblock_size,
			 q->unblock_buffer, 0);
  q->req_index = 0;
  q->req_index_addr = (uint64_t) & q->req_index;

  q->req_processed = 0;
  q->req_processed_reset_limit = q->req_count;
  bgpdcmf_mem_hdl->type = REQUEST_BUFFER;
  DCMF_Memregion_create (&q->req_index_mem_hdl, &bytes,
			 sizeof (uint64_t), (void *) q->req_index_addr, 0);
  for (i = 0; i < MAX_CONNECTION; ++i)
    q->unblock_buffer[i].req_rank = -1;

  start_connection_listener_thread ();
  start_index_thread ();
  return 0;
}

static int
server_req_queue_destroy (bgpdcmf_request_queue_handle * q)
{
  return 0;
}

void
dump_packet (char *buf)
{

  char *
    pkt = buf;
  int
    i;
  for (i = 0; i < 1200; i++)
    fprintf (stderr, "%c", pkt[i]);
  fprintf (stderr, "\n");
}

static NNTI_result_t
insert_buf_bufhash (NNTI_buffer_t * buf)
{
  NNTI_result_t
    rc = NNTI_OK;
  uint32_t
    h =
    hash6432shift ((uint64_t) buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.
		   buf);

  nthread_lock (&nnti_buf_bufhash_lock);
  assert (buffers_by_bufhash.find (h) == buffers_by_bufhash.end ());
  buffers_by_bufhash[h] = buf;
  nthread_unlock (&nnti_buf_bufhash_lock);

  //log_debug (nnti_debug_level, "bufhash buffer added (buf=%p)", buf);

  return (rc);
}

static NNTI_buffer_t *
get_buf_bufhash (const uint32_t bufhash)
{
  NNTI_buffer_t *
    buf = NULL;
/*
  log_debug (nnti_debug_level, "looking for bufhash=%llu",
	     (uint64_t) bufhash);
*/
  nthread_lock (&nnti_buf_bufhash_lock);
  if (buffers_by_bufhash.find (bufhash) != buffers_by_bufhash.end ())
    {
      buf = buffers_by_bufhash[bufhash];
    }
  nthread_unlock (&nnti_buf_bufhash_lock);

  if (buf != NULL)
    {
      //log_debug (nnti_debug_level, "buffer found (buf=%p)", buf);
      return buf;
    }
  return (NULL);
}

static NNTI_buffer_t *
del_buf_bufhash (NNTI_buffer_t * buf)
{
  uint32_t
    h =
    hash6432shift ((uint64_t) buf->buffer_addr.NNTI_remote_addr_t_u.bgpdcmf.
		   buf);
  log_level
    debug_level = nnti_debug_level;

  nthread_lock (&nnti_buf_bufhash_lock);
  if (buffers_by_bufhash.find (h) != buffers_by_bufhash.end ())
    {
      buf = buffers_by_bufhash[h];
    }
  nthread_unlock (&nnti_buf_bufhash_lock);

  if (buf != NULL)
    {
      log_debug (debug_level, "buffer found");
      buffers_by_bufhash.erase (h);
    }
  else
    {
      log_debug (debug_level, "buffer NOT found");
    }

  return (buf);
}

static void
print_bufhash_map ()
{
  if (!logging_debug (nnti_debug_level))
    {
      return;
    }

  buf_by_bufhash_iter_t
    i;
  for (i = buffers_by_bufhash.begin (); i != buffers_by_bufhash.end (); i++)
    {
      log_debug (nnti_debug_level, "bufhash_map key=%llu buf=%p", i->first,
		 i->second);
    }
}

static NNTI_result_t
insert_wr_wrhash (dcmf_work_request * wr)
{
  NNTI_result_t
    rc = NNTI_OK;
  uint32_t
    h = hash6432shift ((uint64_t) wr);

  nthread_lock (&nnti_wr_wrhash_lock);
  assert (wr_by_wrhash.find (h) == wr_by_wrhash.end ());
  wr_by_wrhash[h] = wr;
  nthread_unlock (&nnti_wr_wrhash_lock);

  log_debug (nnti_debug_level,
	     "wrhash work request added (wr=%p ; wr.hash=%llu)", wr,
	     (uint64_t) h);

  return (rc);
}

static dcmf_work_request *
get_wr_wrhash (const uint32_t wrhash)
{
  dcmf_work_request *
    wr = NULL;

  //log_debug (nnti_debug_level, "looking for wrhash=%llu", (uint64_t) wrhash);
  nthread_lock (&nnti_wr_wrhash_lock);
  if (wr_by_wrhash.find (wrhash) != wr_by_wrhash.end ())
    {
      wr = wr_by_wrhash[wrhash];
    }
  nthread_unlock (&nnti_wr_wrhash_lock);

  if (wr != NULL)
    {
      log_debug (nnti_debug_level, "work request found (wr=%p)", wr);
      return wr;
    }

 log_debug (nnti_debug_level, "work request NOT found");
  print_wrhash_map ();

  return (NULL);
}

static dcmf_work_request *
del_wr_wrhash (dcmf_work_request * wr)
{
  uint32_t
    h = hash6432shift ((uint64_t) wr);
  log_level
    debug_level = nnti_debug_level;

  nthread_lock (&nnti_wr_wrhash_lock);
  if (wr_by_wrhash.find (h) != wr_by_wrhash.end ())
    {
      wr = wr_by_wrhash[h];
    }
  nthread_unlock (&nnti_wr_wrhash_lock);

  if (wr != NULL)
    {
      //log_debug (debug_level, "work request found wr=%p", wr);
      wr_by_wrhash.erase (h);
    }
  else
    {
      log_debug (debug_level, "work request NOT found");
    }

  return (wr);
}

static void
print_wrhash_map ()
{
  if (!logging_debug (nnti_debug_level))
    {
      return;
    }

  wr_by_wrhash_iter_t
    i;
  for (i = wr_by_wrhash.begin (); i != wr_by_wrhash.end (); i++)
    {
      log_debug (nnti_debug_level, "wrhash_map key=%llu wr=%p", i->first,
		 i->second);
    }
}


static int8_t
is_any_buf_op_complete (const NNTI_buffer_t ** buf_list,
			const uint32_t buf_count, uint32_t * which)
{
  int8_t
    rc = FALSE;

  log_debug (nnti_debug_level, "enter");

  for (uint32_t i = 0; i < buf_count; i++)
    {
      if ((buf_list[i] != NULL) &&
	  (is_wr_queue_empty (buf_list[i]) == FALSE) &&
	  (is_buf_op_complete (buf_list[i]) == TRUE))
	{

	  *which = i;
	  rc = TRUE;
	  break;
	}
    }

  log_debug (nnti_debug_level, "exit (rc=%d)", rc);

  return (rc);
}

static int8_t
is_wr_queue_empty (const NNTI_buffer_t * reg_buf)
{
  int8_t
    rc = FALSE;
  bgpdcmf_memory_handle *
    dcmf_mem_hdl = NULL;

  log_debug (nnti_debug_level, "enter");

  dcmf_mem_hdl = (bgpdcmf_memory_handle *) reg_buf->transport_private;
  assert (dcmf_mem_hdl);

  if (dcmf_mem_hdl->wr_queue.empty ())
    {
      rc = TRUE;
    }

  log_debug (nnti_debug_level, "exit (rc=%d)", rc);
  return (rc);
}


static int8_t
is_buf_op_complete (const NNTI_buffer_t * reg_buf)
{
  int8_t
    rc = FALSE;
  bgpdcmf_memory_handle *
    dcmf_mem_hdl = NULL;
  dcmf_work_request *
    wr = NULL;

  log_debug (nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

  dcmf_mem_hdl = (bgpdcmf_memory_handle *) reg_buf->transport_private;

  if (is_wr_queue_empty (reg_buf) == TRUE)
    {
      log_debug (nnti_debug_level,
		 "work request queue is empty - return FALSE");
      rc = FALSE;
    }
  else
    {
      wr = dcmf_mem_hdl->wr_queue.front ();
      assert (wr);

      rc = is_wr_complete (wr);
    }

  if (rc == TRUE)
    {
      log_debug (nnti_debug_level, "op is complete");
    }
  log_debug (nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

  return (rc);
}

static int8_t
is_wr_complete (dcmf_work_request * wr)
{
  int8_t
    rc = FALSE;
  bgpdcmf_memory_handle *
    dcmf_mem_hdl = NULL;

  dcmf_mem_hdl = (bgpdcmf_memory_handle *) wr->reg_buf->transport_private;
  assert (dcmf_mem_hdl);

  switch (dcmf_mem_hdl->type)
    {
    case SEND_BUFFER:
      if (wr->op_state == SEND_COMPLETE)
	{
	  rc = TRUE;
	}
      break;
    case PUT_SRC_BUFFER:
      if (wr->op_state == RDMA_COMPLETE)
	{
	  rc = TRUE;
	}
      break;
    case GET_DST_BUFFER:
      if (wr->op_state == RDMA_COMPLETE)
	{
	  rc = TRUE;
	}
      break;
    case REQUEST_BUFFER:
    case RECEIVE_BUFFER:
      if (wr->op_state == RECV_COMPLETE)
	{
	  rc = TRUE;
	}
      break;
    case PUT_DST_BUFFER:
      if (wr->op_state == RDMA_COMPLETE)
	{
	  rc = TRUE;
	}
      break;
    case GET_SRC_BUFFER:
      if (wr->op_state == RDMA_COMPLETE)
	{
	  rc = TRUE;
	}
      break;
    case RDMA_TARGET_BUFFER:
      if (wr->op_state == RDMA_COMPLETE)
	{
	  rc = TRUE;
	}
      break;
    case UNKNOWN_BUFFER:
      break;
    }

  log_debug (nnti_debug_level, "exit (rc=%d)", rc);
  return (rc);
}

static NNTI_result_t
repost_recv_work_request (NNTI_buffer_t * reg_buf, dcmf_work_request * wr)
{
  bgpdcmf_memory_handle *
    bgpdcmf_mem_hdl = NULL;

  log_debug (nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

  bgpdcmf_mem_hdl = (bgpdcmf_memory_handle *) reg_buf->transport_private;
  assert (bgpdcmf_mem_hdl);

  //memset (&wr->wc, 0, sizeof (nnti_bgpdcmf_work_completion));

  wr->last_op = 0;
  wr->wc.op_complete = 0 ;

  //memset (&wr->op_state, 0, sizeof (bgpdcmf_op_state_t));
  //register_wc (wr);

  bgpdcmf_mem_hdl->wr_queue.push_back(wr);

  log_debug (nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

  return (NNTI_OK);
}


