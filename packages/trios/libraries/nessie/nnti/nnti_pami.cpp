/*
 * nnti_pami.c
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
#include <stdint.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <pami.h>
#include <hwi/include/bqc/A2_inlines.h>

#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <assert.h>

#include <sys/time.h>
#include <sys/resource.h>

#include <map>
#include <deque>

#include "nnti_pami.h"
#include "nnti_utils.h"
#include "nnti_internal.h"

#define TRACE_ERR(x)		/* fprintf  x */

/* GLOBAL VAR */


#define BUFSIZE  2048
#define ITERATIONS  1000
#define MAX_CONNECTION  1024
size_t my_rank;
int is_client = 0;
volatile unsigned _recv_active;
volatile unsigned send_active;
volatile unsigned _recv_iteration;
char _recv_buffer[BUFSIZE] __attribute__ ((__aligned__ (16)));
volatile static int insert_conn_index;
int conn_list[1024][2];  /*Populate this with rank and need Attr during dispatch_recv immediate */
int cur_conn_idx =0;
int conn_ack =0;
pami_context_t * contexts;
size_t dispatch_id                 = 37;
pthread_t Progress_thread;
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
} bgqpami_connection_state;

typedef enum
{
  SERVER_CONNECTION,
  CLIENT_CONNECTION
} bgqpami_connection_type;

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
} bgqpami_buffer_type;


#define PAMI_OP_PUT_INITIATOR  1
#define PAMI_OP_GET_INITIATOR  2
#define PAMI_OP_PUT_TARGET     3
#define PAMI_OP_GET_TARGET     4
#define PAMI_OP_SEND_REQ           5
#define PAMI_OP_NEW_REQUEST    6
#define PAMI_OP_RESULT         7
#define PAMI_OP_RECEIVE        8
#define PAMI_OP_REGISTER_RDMA        9

typedef enum
{
  SEND_COMPLETE,
  RECV_COMPLETE,
  RDMA_WRITE_INIT,
  RDMA_TARGET_INIT,
  RDMA_READ_INIT,
  RDMA_COMPLETE
} bgqpami_op_state_t;


typedef struct
{

  size_t req_rank;
  pami_memregion_t response_hdl;
} nnti_bgqpami_fetch_index_req;

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
  uint64_t  dest_addr;
  uint16_t op;
} nnti_bgqpami_work_completion;

/**
 * attrs to send to the client
 */
typedef struct
{
  uint64_t req_index_addr;	/* address of the request index var on the server */
  pami_memregion_t req_index_mem_hdl;

  uint64_t req_buffer_addr;	/* address of the request buffer */
  uint64_t req_size;		/* the maximum size of a request */
  uint64_t req_count;		/* the number of requests that will fit in the queue */
  pami_memregion_t req_mem_hdl;
  uint64_t wc_buffer_addr;	/* address of the work completion array */
  pami_memregion_t wc_mem_hdl;
  uint64_t unblock_buffer_addr;	/* address of the fetch_req_array  */
  pami_memregion_t unblock_mem_hdl;
} nnti_bgqpami_server_queue_attrs;

/**
 * attrs to send to the server
 */
typedef struct
{
  uint64_t req_index;		/* this buffer contains the current index of the request queue */
  uint64_t req_index_addr;	/* address of the index buffer */
  pami_memregion_t req_index_mem_hdl;	/* mem hdl to send to the server */
  pami_memregion_t fetch_index_hdl;
} nnti_bgqpami_client_queue_attrs;

typedef union
{
  nnti_bgqpami_client_queue_attrs client;
  nnti_bgqpami_server_queue_attrs server;
} nnti_bgqpami_queue_remote_attrs;

typedef struct
{
  uint32_t peer_rank;
  nnti_bgqpami_client_queue_attrs queue_local_attrs;
  nnti_bgqpami_queue_remote_attrs queue_remote_attrs;
  bgqpami_connection_state state;
  bgqpami_connection_type connection_type;
} bgqpami_connection;


typedef struct
{
  uint8_t is_initiator;

  const NNTI_buffer_t *reg_buf;
  nnti_bgqpami_work_completion wc;
  uint64_t wc_dest_addr;
  pami_memregion_t wc_dest_mem_hdl;
  pami_memregion_t wc_mem_hdl;
  uint8_t last_op;
  bgqpami_op_state_t op_state;
  uint64_t is_last_op_complete;
  uint32_t peer_rank;
} pami_work_request;


typedef
  std::deque <
pami_work_request * >
  wr_queue_t;
typedef
  std::deque <
pami_work_request * >::iterator
  wr_queue_iter_t;

typedef struct
{
  bgqpami_buffer_type
    type;
  pami_memregion_t
    mem_hdl;			/* actual data memory handle */
  wr_queue_t
    wr_queue;
  nthread_lock_t
    wr_queue_lock;
  uint32_t
    ref_count;
} bgqpami_memory_handle;

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
  pami_memregion_t
    req_index_mem_hdl;

  char *
    req_buffer;			/* pointer to the head of the request buffer */
  uint64_t
    req_size;			/* the maximum size of a request */
  uint64_t
    req_count;			/* the number of requests that will fit in the queue */
  uint64_t
    req_buffer_size;		/* the size of the request buffer in bytes (req_size*req_count) */
  pami_memregion_t
    req_mem_hdl;

  nnti_bgqpami_work_completion *
    wc_buffer;			/* pointer to work completion array */
  uint64_t
    wc_buffer_size;		/* the size of the work completion buffer in bytes  */
  pami_memregion_t
    wc_mem_hdl;
  nnti_bgqpami_fetch_index_req *
    unblock_buffer;		/* pointer to individual fetch index array */
  pami_memregion_t
    unblock_mem_hdl;
  uint64_t
    req_processed_reset_limit;
  uint64_t
    req_processed;
  uint64_t
    total_req_processed;
} bgqpami_request_queue_handle;

typedef struct
{
  uint32_t
    myrank;
  uint32_t
    remote_rank;
  int
    mypid;
  bgqpami_request_queue_handle
    req_queue;
} bgqpami_transport_global;

static nthread_lock_t
  nnti_bgqpami_lock;
static nthread_lock_t
  nnti_index_lock;
static int
start_connection_listener_thread (void);
static NNTI_result_t
register_memory (bgqpami_memory_handle * hdl,
		 void *buf, uint64_t len, const NNTI_buf_ops_t remote_op);
static int
unregister_memory (bgqpami_memory_handle * hdl);
static NNTI_result_t
process_event (const NNTI_buffer_t * reg_buf,
	       const NNTI_buf_ops_t remote_op,
	       pami_work_request * wr, const int timeout);
/* NEW SHYAMALI */

static NNTI_result_t
repost_recv_work_request (NNTI_buffer_t * reg_buf, pami_work_request * wr);

static pami_work_request *
first_incomplete_wr (bgqpami_memory_handle * _hdl);
static int8_t
is_buf_op_complete (const NNTI_buffer_t * reg_buf);
static int8_t
is_wr_complete (pami_work_request * wr);
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
recv_full (int rank, void *buf, size_t num);
static NNTI_result_t
init_connection (bgqpami_connection ** conn, int rank, const int is_server);
static void
close_connection (bgqpami_connection * c);

static NNTI_result_t
insert_conn_rank (const uint32_t rank, bgqpami_connection * conn);
static bgqpami_connection *
get_conn_rank (const uint32_t rank);
static bgqpami_connection *
del_conn_rank (const uint32_t rank);

static int
server_req_queue_init (bgqpami_request_queue_handle * q,
		       char *buffer, uint64_t req_size, uint64_t req_count);
static int
server_req_queue_destroy (bgqpami_request_queue_handle * q);

static int
client_req_queue_init (bgqpami_connection * c);
static int
client_req_queue_destroy (bgqpami_connection * c);

static int
reset_req_index (bgqpami_request_queue_handle * req_queue_attrs);

static int
fetch_server_req_buffer_offset (nnti_bgqpami_client_queue_attrs *
				local_req_queue_attrs,
				nnti_bgqpami_server_queue_attrs *
				remote_req_queue_attrs,
				uint64_t addend, uint64_t * offset, int rank);
static int
send_req (nnti_bgqpami_client_queue_attrs * local_req_queue_attrs,
	  nnti_bgqpami_server_queue_attrs * remote_req_queue_attrs,
	  uint64_t offset, const NNTI_buffer_t * reg_buf,
	  int rank, pami_work_request * wr);

static  int
send_req_wc (nnti_bgqpami_client_queue_attrs * local_req_queue_attrs,
	     nnti_bgqpami_server_queue_attrs * remote_req_queue_attrs,
	     const NNTI_buffer_t * reg_buf, uint64_t offset, pami_work_request * wr, int rank);

static int
request_send (nnti_bgqpami_client_queue_attrs * client_q,
	      nnti_bgqpami_server_queue_attrs * server_q,
	      const NNTI_buffer_t * reg_buf, int rank);

/* NEW */
static NNTI_result_t
register_wc (pami_work_request * wr);
static NNTI_result_t
insert_buf_bufhash (NNTI_buffer_t * buf);
static NNTI_buffer_t *
get_buf_bufhash (const uint32_t bufhash);
static NNTI_buffer_t *
del_buf_bufhash (NNTI_buffer_t * buf);
static void
print_bufhash_map (void);

static NNTI_result_t
insert_wr_wrhash (pami_work_request *);
static pami_work_request *
get_wr_wrhash (const uint32_t bufhash);
static pami_work_request *
del_wr_wrhash (pami_work_request *);
static void
print_wrhash_map (void);
static void
send_rdma_wc (pami_work_request * wr,
	      const NNTI_buffer_t * local_buf,
	      const NNTI_buffer_t * remote_buf);


static
  NNTI_result_t
new_server_connection (bgqpami_connection * c, int peer_rank);
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
bgqpami_connection * >connections_by_rank;
typedef
  std::map < int,
bgqpami_connection * >::iterator
  conn_by_rank_iter_t;
typedef
  std::pair < int,
bgqpami_connection * >
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
pami_work_request * >
  wr_by_wrhash;
typedef
  std::map <
  uint32_t,
pami_work_request * >::iterator
  wr_by_wrhash_iter_t;
typedef
  std::pair <
  uint32_t,
pami_work_request * >
  wr_by_wrhash_t;
static nthread_lock_t
  nnti_wr_wrhash_lock;



/* NEW */

static bgqpami_transport_global
  transport_global_data;
static const int
  MIN_TIMEOUT = 80000;		/* in milliseconds */
static log_level
  nnti_event_debug_level;
static nnti_bgqpami_fetch_index_req
  client_fetch_req;

/**************************  Callback functions for PAMI_Send / Recv protocol ************/
/* Low overhead timing function */

static void * Progress_function(void * dummy)
{
        pami_result_t result = PAMI_ERROR;

        while (1)
        {
        result = PAMI_Context_trylock_advancev(&(contexts[1]), 5, 1000);
                usleep(100);
        }

        return NULL;
}



void *
safemalloc(size_t size)
{
  void *
    ptr;
  int
    rc;

  rc = posix_memalign (&ptr, 128, size);
  if (rc)
    return NULL;

  return ptr;
}


void dispatch_done_cb(pami_context_t context, void * cookie, pami_result_t result)
{

  fflush(stdout);
  return;
}
static void cb_done (void *ctxt, void * clientdata, pami_result_t err)
{
  int * active = (int *) clientdata;
  (*active)--;
}

void dispatch_recv_cb(pami_context_t context,
                      void * cookie,
                      const void * header_addr, size_t header_size,
                      const void * pipe_addr,
                      size_t data_size,
                      pami_endpoint_t origin,
                      pami_recv_t * recv)
{
  pami_result_t result = PAMI_ERROR;
  bgqpami_connection *conn;

  size_t task;
  size_t ctxoff;
  result = PAMI_Endpoint_query(origin, &task, &ctxoff);

   if ((data_size < 32 ) || (data_size > 80 ))
	return ;
  conn_ack--;

  if (pipe_addr!=NULL)
  {

    	nthread_lock (&nnti_bgqpami_lock);
        conn = get_conn_rank (origin);
  	if (conn == NULL) {
        	if(is_client != 1){
			conn = (bgqpami_connection *) calloc (1, sizeof (bgqpami_connection));
          conn->peer_rank = origin;
	  memcpy((void *)&conn->queue_remote_attrs.client, pipe_addr, data_size);
          insert_conn_rank (conn->peer_rank, conn);

          log_debug (nnti_debug_level,
                     "accepted new connection from 192.168.1.%d", conn->peer_rank);
		   conn_list[cur_conn_idx][0] = origin;
		   conn_list[cur_conn_idx][1] = 0;
		   cur_conn_idx++;
		}
  	}
	else if((is_client == 1) &&  (conn != NULL)) {
	  memcpy((void *)&conn->queue_remote_attrs.server, pipe_addr, data_size);
	  conn->peer_rank = origin;
	  
         }	
  	nthread_unlock (&nnti_bgqpami_lock);
  }
  else
  {
    fflush(stdout);
    nthread_lock (&nnti_bgqpami_lock);
    conn = get_conn_rank (origin);
     if (conn == NULL) {
          conn = (bgqpami_connection *) calloc (1, sizeof (bgqpami_connection));
      }
          conn->peer_rank = origin;
          insert_conn_rank (conn->peer_rank, conn);

    recv->cookie      = (void *) conn;
    recv->local_fn    = dispatch_done_cb;
    if(is_client == 1)
    		recv->addr        = &conn->queue_remote_attrs.server;
    recv->type        = PAMI_TYPE_BYTE;
    recv->offset      = 0;
    recv->data_fn     = PAMI_DATA_COPY;
    recv->data_cookie = NULL;
    nthread_unlock (&nnti_bgqpami_lock);
  }

  return;
}

static void decrement(void *ctxt, void * clientdata, pami_result_t err)
{
  int * active = (int *) clientdata;
  (*active)--;
}


void
send_once (char *buffer, size_t sndlen, size_t targetrank)
{
  pami_endpoint_t target_ep;
   pami_result_t  result = PAMI_Endpoint_create(client, (pami_task_t) targetrank, 1, &target_ep);
  
                                     // void *iov_base - Contains the address of a buffer. 
                                     // size_t iov_len - Contains the length of the buffer.

                                        int header  = 37373;

                                           pami_send_t parameters;
                                          parameters.send.header.iov_base = &header;
                                          parameters.send.header.iov_len  = sizeof(int);
                                          parameters.send.data.iov_base   = buffer;
                                          parameters.send.data.iov_len    = sndlen;
                                          parameters.send.dispatch        = dispatch_id;
                                          parameters.send.dest            = target_ep;
                                          parameters.events.cookie        = (void *)&_recv_active;
                                          parameters.events.local_fn      = cb_done;
                                          parameters.events.remote_fn     = NULL;//cb_done;
  
                                          uint64_t t0 = GetTimeBase();
			                  _recv_active = 1;	

                                          result = PAMI_Send(contexts[0], &parameters);

                                          while (_recv_active)
                                          {
                                           result = PAMI_Context_trylock_advancev(&(contexts[0]), 1, 1000);
                                          }
}



/* ==================== */




/**
 * @brief Initialize NNTI to use a specific transport.
 *
 */
NNTI_result_t
NNTI_bgqpami_init (const NNTI_transport_id_t trans_id,
		   const char *my_url, NNTI_transport_t * trans_hdl)
{
  static int
    initialized = 0;

  NNTI_result_t
    rc = NNTI_OK;

  pami_result_t result = PAMI_ERROR;
  size_t world_size, world_rank;
  int i;

  /* initialize the second client */
  char  clientname[20] = "NNTI";
  pami_client_t client;
  result = PAMI_Client_create((char *)&clientname, &client, NULL, 0);

  /* query properties of the client */
  pami_configuration_t config[3];
  size_t num_contexts;

  config[0].name = PAMI_CLIENT_NUM_TASKS;
  config[1].name = PAMI_CLIENT_TASK_ID;
  config[2].name = PAMI_CLIENT_NUM_CONTEXTS;
  result = PAMI_Client_query(client, config, 3);
  world_size   = config[0].value.intval;
  my_rank   = config[1].value.intval;
  num_contexts = config[2].value.intval;


  /* initialize the contexts */
  contexts = (pami_context_t *) safemalloc( num_contexts * sizeof(pami_context_t) );

  result = PAMI_Context_createv(client, NULL, 0, contexts, num_contexts );
  if (my_rank==0)
  {
    printf("hello world from rank %ld of %ld no. of contexts %ld\n", my_rank, world_size, num_contexts);
    fflush(stdout);
  }



  initialized = 0;
  _recv_active = 1;
  send_active = 1;
  _recv_iteration = 0;

  if (!initialized)
    {
      for (i=0; i< 1024; i++)
	{
	 conn_list[i][0]= -1;
	 conn_list[i][1]= -1;
	}
      nthread_lock_init (&nnti_bgqpami_lock);
      nthread_lock_init (&nnti_index_lock);
      nthread_lock_init (&nnti_conn_peer_lock);
      nthread_lock_init (&nnti_wr_wrhash_lock);
      nthread_lock_init (&nnti_buf_bufhash_lock);
      pami_dispatch_callback_function dispatch_cb;
  dispatch_cb.p2p                    = dispatch_recv_cb;
  pami_dispatch_hint_t dispatch_hint = {0};
   _recv_active = 1;
  /* setting this hint to "disable" will disable the pipe_addr!=NULL code path */
  //dispatch_hint.recv_immediate       = PAMI_HINT_DISABLE; 
    result = PAMI_Dispatch_set(contexts[0], dispatch_id, dispatch_cb, (void *)&_recv_active, dispatch_hint);
      result = PAMI_Dispatch_set(contexts[1], dispatch_id, dispatch_cb, (void *)&_recv_active, dispatch_hint);

      memset (&transport_global_data, 0, sizeof (bgqpami_transport_global));

      log_debug (nnti_debug_level, "my_url=%s", my_url);


      log_debug (nnti_debug_level, "initializing Blue Gene  DMA PAMI device");

      transport_global_data.mypid = getpid ();
      transport_global_data.myrank = my_rank;


      create_peer (&trans_hdl->me, transport_global_data.myrank);
      int status = pthread_create(&Progress_thread, NULL, &Progress_function, NULL);
      assert(status==0);

      initialized = TRUE;
    }

  log_debug (nnti_debug_level, "exit");

  return (rc);
}


/**
 * @brief Return the URL field of this transport.
 *
 * Return the URL field of this transport.  After initialization, the transport will
 * have a specific location on the network where peers can contact it.  The
 * transport will convert this location to a string that other instances of the
 * transport will recobgqpamize.
 *
 */
NNTI_result_t
NNTI_bgqpami_get_url (const NNTI_transport_t * trans_hdl,
		      char *url, const uint64_t maxlen)
{
  NNTI_result_t
    rc = NNTI_OK;

  assert (trans_hdl);
  assert (url);
  assert (maxlen > 0);

  log_debug (nnti_debug_level, "enter");

  strncpy (url, trans_hdl->me.url, maxlen);
  url[maxlen - 1] = '\0';

  log_debug (nnti_debug_level, "exit");

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
NNTI_bgqpami_connect (const NNTI_transport_t * trans_hdl,
		      const char *url,
		      const int timeout, NNTI_peer_t * peer_hdl)
{
  NNTI_result_t
    rc = NNTI_OK;

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

  bgqpami_connection *
    c = NULL;


  assert (trans_hdl);
  assert (peer_hdl);

  log_debug (nnti_debug_level, "enter");
  is_client = 1;

  if (url != NULL)
    {
      if ((rc =
	   nnti_url_get_transport (url, transport, NNTI_URL_LEN)) != NNTI_OK)
	{
	  return (rc);
	}
      if (0 != strcmp (transport, "pami"))
	{
	  /* the peer described by 'url' is not a BGP Torus peer */
	  return (NNTI_EINVAL);
	}

      if ((rc = nnti_url_get_address (url, address, NNTI_URL_LEN)) != NNTI_OK)
	{
	  return (rc);
	}
      temp = strdup (url);
      tmp = strtok (temp, "//");
      pset_rank = strtol (tmp + 1, NULL, 0);
      transport_global_data.remote_rank = pset_rank;	/* This is server rank */
    }
  else
    {
      /*  */
      rc = NNTI_EINVAL;
      goto cleanup;
    }
    nnti_bgqpami_client_queue_attrs  ca_out;

  c = (bgqpami_connection *) calloc (1, sizeof (bgqpami_connection));
  conn_ack = 1;

  client_req_queue_init (c);
  memset (&ca_out, 0, sizeof (ca_out));
  ca_out.req_index = c->queue_local_attrs.req_index;
  ca_out.req_index_addr = c->queue_local_attrs.req_index_addr;
  memcpy (&ca_out.req_index_mem_hdl,
          &c->queue_local_attrs.req_index_mem_hdl, sizeof (pami_memregion_t));
  send_once ((char *)&ca_out, sizeof(ca_out), transport_global_data.remote_rank);
  insert_conn_rank (pset_rank, c);

  if (c == NULL)
    {
      rc = NNTI_EIO;
      goto cleanup;
    }
   while (conn_ack > 0){
   	PAMI_Context_trylock_advancev(&(contexts[1]), 1, 1000);
    }  
  create_peer (peer_hdl, pset_rank);
cleanup:
  log_debug (nnti_debug_level, "exit");
  return (rc);
}


/**
 * @brief Terminate communication with this peer.
 *
 * Perform any transport specific actions necessary to end communication with
 * this peer.
 */
NNTI_result_t
NNTI_bgqpami_disconnect (const NNTI_transport_t * trans_hdl,
			 NNTI_peer_t * peer_hdl)
{
  NNTI_result_t
    rc = NNTI_OK;

  assert (trans_hdl);
  assert (peer_hdl);

  log_debug (nnti_debug_level, "enter");
  bgqpami_connection *
    conn =
    get_conn_rank (peer_hdl->peer.NNTI_remote_process_t_u.bgqpami.pset_rank);
  close_connection (conn);
  del_conn_rank (peer_hdl->peer.NNTI_remote_process_t_u.bgqpami.pset_rank);


  log_debug (nnti_debug_level, "exit");

  return (rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 */
NNTI_result_t
NNTI_bgqpami_register_memory (const NNTI_transport_t * trans_hdl,
			      char *buffer,
			      const uint64_t size,
			      const uint64_t num_elements,
			      const NNTI_buf_ops_t ops,
			      const NNTI_peer_t * peer,
			      NNTI_buffer_t * reg_buf)
{
  NNTI_result_t
    rc = NNTI_OK;

  pami_work_request *
    wr = NULL;
  NNTI_buffer_t *
    old_buf = NULL;
  bgqpami_memory_handle *
    bgqpami_mem_hdl = NULL;


  assert (trans_hdl);
  assert (buffer);
  assert (size > 0);
  assert (ops > 0);
  assert (reg_buf);

  log_debug (nnti_debug_level, "enter");

  old_buf = get_buf_bufhash (hash6432shift ((uint64_t) buffer));
  if (old_buf == NULL)
    {
      bgqpami_mem_hdl = new bgqpami_memory_handle ();
      bgqpami_mem_hdl->ref_count = 1;
      nthread_lock_init (&bgqpami_mem_hdl->wr_queue_lock);
    }
  else
    {
      bgqpami_mem_hdl = (bgqpami_memory_handle *) old_buf->transport_private;
      bgqpami_mem_hdl->ref_count++;
    }

  assert (bgqpami_mem_hdl);

  reg_buf->transport_id = trans_hdl->id;
  reg_buf->buffer_owner = trans_hdl->me;
  reg_buf->ops = ops;
  reg_buf->payload_size = size;
  reg_buf->payload = (uint64_t) buffer;
  reg_buf->transport_private = (uint64_t) bgqpami_mem_hdl;

  if (peer != NULL)
    {
      reg_buf->buffer_owner = *peer;
    }


  log_debug (nnti_debug_level, "rpc_buffer->payload_size=%ld",
	     reg_buf->payload_size);
  reg_buf->buffer_addr.transport_id = NNTI_TRANSPORT_PAMI;
  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.size = size;
  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.buf = (uint64_t) buffer;
  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.owner_rank =
    my_rank;
  if (bgqpami_mem_hdl->ref_count == 1)
    {
      if (ops == NNTI_RECV_QUEUE)
	{
	  bgqpami_request_queue_handle *
	    q_hdl = &transport_global_data.req_queue;

	  /*
	   * This is a receive-only buffer.  This buffer is divisible by
	   * NNTI_REQUEST_BUFFER_SIZE.  This buffer can hold more than
	   * one short request.  Assume this buffer is a request queue.
	   */
	  bgqpami_mem_hdl->type = REQUEST_BUFFER;
	  memset (q_hdl, 0, sizeof (bgqpami_request_queue_handle));
	  q_hdl->reg_buf = reg_buf;	/* req buffer can be accessed from global req_queue  */
	  server_req_queue_init (q_hdl, buffer, size, num_elements);

	  reg_buf->payload_size = q_hdl->req_size;
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.type =
	    NNTI_PAMI_REQUEST_BUFFER;

	}
      else if (ops == NNTI_RECV_DST)
	{
	  if (size == NNTI_RESULT_BUFFER_SIZE)
	    {
	      /*
	       * This is a receive-only buffer.  This buffer can hold exactly
	       * one short result.  Assume this buffer is a result queue.
	       */
	      bgqpami_mem_hdl->type = RESULT_BUFFER;

	      rc =
		register_memory (bgqpami_mem_hdl, buffer,
				 NNTI_RESULT_BUFFER_SIZE, ops);

	      reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.type =
		NNTI_PAMI_RECEIVE_DST;
	    }
	  else
	    {
	      /*
	       * This is a receive-only buffer.  This buffer doesn't look
	       * like a request buffer or a result buffer.  I don't know
	       * what it is.  Assume it is a regular data buffer.
	       */
	      bgqpami_mem_hdl->type = RECEIVE_BUFFER;

	      rc = register_memory (bgqpami_mem_hdl, buffer, size, ops);
	      reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.type =
		NNTI_PAMI_RECEIVE_DST;

	    }

	}
      else if (ops == NNTI_SEND_SRC)
	{
	  bgqpami_mem_hdl->type = SEND_BUFFER;

	  rc = register_memory (bgqpami_mem_hdl, buffer, size, ops);
	  if (rc != NNTI_OK)
	    {
	      fprintf (stderr,
		       "failed registering short request in register memory\n");
	    }

	}
      else if (ops == NNTI_GET_DST)
	{
	  bgqpami_mem_hdl->type = GET_DST_BUFFER;
	  rc = register_memory (bgqpami_mem_hdl, buffer, size, ops);

	}
      else if (ops == NNTI_GET_SRC)
	{
	  bgqpami_mem_hdl->type = GET_SRC_BUFFER;


	  rc = register_memory (bgqpami_mem_hdl, buffer, size, ops);

	}
      else if (ops == NNTI_PUT_SRC)
	{
	  bgqpami_mem_hdl->type = PUT_SRC_BUFFER;


	  rc = register_memory (bgqpami_mem_hdl, buffer, size, ops);
	}
      else if (ops == NNTI_PUT_DST)
	{
	  bgqpami_mem_hdl->type = PUT_DST_BUFFER;

	  rc = register_memory (bgqpami_mem_hdl, buffer, size, ops);

	}
      else if (ops == (NNTI_GET_SRC | NNTI_PUT_DST))
	{
	  bgqpami_mem_hdl->type = RDMA_TARGET_BUFFER;

	  rc = register_memory (bgqpami_mem_hdl, buffer, size, ops);
	}
      else
	{
	  bgqpami_mem_hdl->type = UNKNOWN_BUFFER;
	}
    }
  wr = (pami_work_request *) calloc (1, sizeof (pami_work_request));
  wr->reg_buf = reg_buf;
  wr->last_op = PAMI_OP_REGISTER_RDMA;
  register_wc (wr);


  if (rc == NNTI_OK)
    {
      reg_buf->buffer_addr.transport_id = NNTI_TRANSPORT_PAMI;
	reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.mem_hdl = (uint64_t)bgqpami_mem_hdl->mem_hdl;
      if (bgqpami_mem_hdl->type == REQUEST_BUFFER)
	{
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.size =
	    transport_global_data.req_queue.req_size;
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.buf =
	    (uint64_t) transport_global_data.req_queue.req_buffer;
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.type =
	    NNTI_PAMI_REQUEST_BUFFER;
	}
      else
	{
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.size =
	    reg_buf->payload_size;
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.buf =
	    (uint64_t) reg_buf->payload;
	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.type =
	    NNTI_PAMI_SEND_SRC;

	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.wc_addr =
	    (uint64_t) & wr->wc;
 	  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.wc_mem_hdl = (uint64_t) &wr->wc_mem_hdl;
	}

    }

  if (bgqpami_mem_hdl->ref_count == 1)
    {
      insert_buf_bufhash (reg_buf);
      log_debug (nnti_debug_level, "bgqpami_mem_hdl->type==%llu",
		 (uint64_t) bgqpami_mem_hdl->type);
      log_debug (nnti_debug_level, "reg_buf.buf.hash==%llu",
		 (uint64_t) hash6432shift (reg_buf->buffer_addr.
					   NNTI_remote_addr_t_u.bgqpami.buf));
    }

  nthread_lock (&bgqpami_mem_hdl->wr_queue_lock);
  bgqpami_mem_hdl->wr_queue.push_back (wr);
  nthread_unlock (&bgqpami_mem_hdl->wr_queue_lock);
  return (rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t
NNTI_bgqpami_unregister_memory (NNTI_buffer_t * reg_buf)
{
  NNTI_result_t
    rc = NNTI_OK;
  bgqpami_memory_handle *
    bgqpami_mem_hdl = NULL;

  assert (reg_buf);

  log_debug (nnti_debug_level, "enter");

  bgqpami_mem_hdl = (bgqpami_memory_handle *) reg_buf->transport_private;

  assert (bgqpami_mem_hdl);
  bgqpami_mem_hdl->ref_count--;

  if (bgqpami_mem_hdl->ref_count == 0)
    {
      log_debug (nnti_debug_level,
		 "bgqpami_mem_hdl->ref_count is 0.  release all resources.");

      if (bgqpami_mem_hdl->type == REQUEST_BUFFER)
	{
	  server_req_queue_destroy (&transport_global_data.req_queue);

	}
      else
	{
	  unregister_memory (bgqpami_mem_hdl);
	}

      del_buf_bufhash (reg_buf);
      nthread_lock (&bgqpami_mem_hdl->wr_queue_lock);
      while (!bgqpami_mem_hdl->wr_queue.empty ())
	{
	  pami_work_request *
	    wr = bgqpami_mem_hdl->wr_queue.front ();
	  log_debug (nnti_debug_level, "removing pending wr=%p", wr);
	  bgqpami_mem_hdl->wr_queue.pop_front ();
	  del_wr_wrhash (wr);
	}
      nthread_unlock (&bgqpami_mem_hdl->wr_queue_lock);
    }
  nthread_lock_fini (&bgqpami_mem_hdl->wr_queue_lock);
  if (bgqpami_mem_hdl)
    delete
      bgqpami_mem_hdl;

  reg_buf->transport_id = NNTI_TRANSPORT_NULL;
  reg_buf->payload_size = 0;
  reg_buf->payload = 0;
  reg_buf->transport_private = 0;

  log_debug (nnti_debug_level, "exit");

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
NNTI_bgqpami_put (const NNTI_buffer_t * src_buffer_hdl,
		  const uint64_t src_offset,
		  const uint64_t src_length,
		  const NNTI_buffer_t * dest_buffer_hdl,
		  const uint64_t dest_offset)
{
  bgqpami_memory_handle *
    src_buf_mem_hdl = NULL;
  NNTI_result_t
    rc = NNTI_OK;
  pami_work_request *
    wr = NULL;
  pami_endpoint_t target_ep;

  src_buf_mem_hdl =
    (bgqpami_memory_handle *) src_buffer_hdl->transport_private;
  size_t
    srv_rank =
    dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgqpami.owner_rank;
  wr = (pami_work_request *) calloc (1, sizeof (pami_work_request));
  assert (wr);
  wr->wc.local_send_complete = 1;
  wr->wc.remote_put_complete = 2;
  pami_result_t result = PAMI_Endpoint_create(client, (pami_task_t) srv_rank, 1, &target_ep);
  pami_put_simple_t parameters;
  parameters.rma.dest                   = target_ep;
  parameters.rma.bytes                  = src_length;
  parameters.rma.cookie                 = (void *)&wr->wc.remote_put_complete;
  parameters.rma.done_fn                = cb_done;
  parameters.addr.local      = (void *)src_buffer_hdl->payload + src_offset;
  parameters.addr.remote     = (void *)dest_buffer_hdl->buffer_addr.
            NNTI_remote_addr_t_u.bgqpami.buf + dest_offset;
  parameters.put.rdone_fn       = cb_done;
  
                   uint64_t t0 = GetTimeBase();
 
                    result = PAMI_Put(contexts[4], &parameters);

                       while (wr->wc.remote_put_complete)
                          {
                                        result = PAMI_Context_trylock_advancev(&(contexts[4]), 1, 1000);
                          }
  
                                               uint64_t t1 = GetTimeBase();
                                                 uint64_t dt = t1-t0;

  register_wc (wr);
  wr->reg_buf = src_buffer_hdl;
  wr->peer_rank = srv_rank;
  memcpy (&wr->wc_dest_mem_hdl,
	  &dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgqpami.
	  wc_mem_hdl, sizeof (pami_memregion_t));
  wr->wc_dest_addr = dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgqpami.wc_addr;
  
  wr->op_state = RDMA_WRITE_INIT;
  wr->wc.op = PAMI_OP_PUT_TARGET;
  wr->wc.req_rank = transport_global_data.myrank;
  wr->wc.byte_len = src_length;
  wr->wc.src_offset = src_offset;
  wr->wc.dest_offset = dest_offset;
  wr->last_op = PAMI_OP_PUT_INITIATOR;
  wr->wc.op_complete = 1;
  wr->wc.local_send_complete = 2; 
  wr->wc.dest_addr = dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgqpami.wc_addr;
  pami_put_simple_t param;
  param.rma.dest                   = target_ep;
  param.rma.bytes                  = sizeof (nnti_bgqpami_work_completion);
  param.rma.cookie                 = (void *)&wr->wc.local_send_complete;
  param.rma.done_fn                = cb_done;
  param.addr.local      = (void *)&wr->wc;
  param.addr.remote     = (void *)dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgqpami.wc_addr;
  param.put.rdone_fn       = cb_done;

                    result = PAMI_Put(contexts[4], &param);

                       while (wr->wc.local_send_complete)
                          {
                                        result = PAMI_Context_trylock_advancev(&(contexts[4]), 1, 1000);
                          }

  nthread_lock (&src_buf_mem_hdl->wr_queue_lock);
  src_buf_mem_hdl->wr_queue.push_back (wr);
  nthread_unlock (&src_buf_mem_hdl->wr_queue_lock);
  insert_wr_wrhash (wr);
  log_debug (nnti_debug_level, "exit");

  return (rc);
}


/**
 * @brief Send a message to a peer.
 *
 * Send a message (<tt>msg_hdl</tt>) to a peer (<tt>peer_hdl</tt>).  It is expected that the
 * message is small, but the exact maximum size is transport dependent.
 */
NNTI_result_t
NNTI_bgqpami_send (const NNTI_peer_t * peer_hdl,
		   const NNTI_buffer_t * msg_hdl,
		   const NNTI_buffer_t * dest_hdl)
{
  NNTI_result_t
    rc = NNTI_OK;


  bgqpami_memory_handle *
    bgqpami_mem_hdl = NULL;

  log_debug (nnti_debug_level, "enter");

  assert (peer_hdl);
  assert (msg_hdl);

  bgqpami_mem_hdl = (bgqpami_memory_handle *) msg_hdl->transport_private;

  if ((dest_hdl == NULL) || (dest_hdl->ops == NNTI_RECV_QUEUE))
    {
      bgqpami_connection *
	conn =
	get_conn_rank (peer_hdl->peer.NNTI_remote_process_t_u.bgqpami.
		       pset_rank);
      assert (conn);
      request_send (&conn->queue_local_attrs,
		    &conn->queue_remote_attrs.server, msg_hdl,
		    conn->peer_rank);

    }
  else
    {
      rc = NNTI_bgqpami_put (msg_hdl, 0, msg_hdl->payload_size, dest_hdl, 0);
      if (rc != NNTI_OK)
	log_error (nnti_debug_level, "Put() failed: %d", rc);
    }

  log_debug (nnti_debug_level, "sending to (%s)", peer_hdl->url);

  log_debug (nnti_debug_level, "exit");

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
NNTI_bgqpami_get (const NNTI_buffer_t * src_buffer_hdl,
		  const uint64_t src_offset,
		  const uint64_t src_length,
		  const NNTI_buffer_t * dest_buffer_hdl,
		  const uint64_t dest_offset)
{
  NNTI_result_t
    rc = NNTI_OK;
  bgqpami_memory_handle *
    target_buf_mem_hdl;
  int
    srv_rank;
  pami_work_request *
    wr = NULL;


  wr = (pami_work_request *) calloc (1, sizeof (pami_work_request));
  wr->reg_buf = dest_buffer_hdl;
  srv_rank =
    src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgqpami.owner_rank;

  target_buf_mem_hdl =
    (bgqpami_memory_handle *) dest_buffer_hdl->transport_private;
  wr->op_state = RDMA_READ_INIT;
  wr->wc.op = PAMI_OP_GET_TARGET;
  wr->wc.req_rank = srv_rank;
  wr->wc.byte_len = src_length;
  wr->wc.src_offset = src_offset;
  wr->wc.dest_offset = dest_offset;
  wr->last_op = PAMI_OP_GET_INITIATOR;
  wr->wc.local_get_complete = 1;
  wr->peer_rank = srv_rank;

  nthread_lock (&nnti_bgqpami_lock);
  pami_endpoint_t target_ep;
  pami_result_t result = PAMI_Endpoint_create(client, (pami_task_t) srv_rank, 1, &target_ep);

  pami_get_simple_t parameters;
  parameters.rma.dest                   = target_ep;
  //parameters.rma.hints          = ;
    parameters.rma.bytes                  = src_length;
      parameters.rma.cookie                 = (void *)&wr->wc.local_get_complete;
       parameters.rma.done_fn                = cb_done;
       parameters.addr.local = (void *)dest_buffer_hdl->payload + dest_offset;
  parameters.addr.remote     = (void *)src_buffer_hdl->buffer_addr.
            NNTI_remote_addr_t_u.bgqpami.buf;
  
                 uint64_t t0 = GetTimeBase();
 
                  result = PAMI_Get(contexts[5], &parameters);

                      while (wr->wc.local_get_complete)
                        {
                                   result = PAMI_Context_trylock_advancev(&(contexts[5]), 1, 1000);
                        }

                                           uint64_t t1 = GetTimeBase();
                                              uint64_t dt = t1-t0;

  memcpy (&wr->wc_dest_mem_hdl,
	  &src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgqpami.
	  wc_mem_hdl, sizeof (pami_memregion_t));
  wr->wc_dest_addr = src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgqpami.wc_addr;
  wr->wc.req_rank = transport_global_data.myrank;
  wr->wc.byte_len = src_length;
  wr->wc.src_offset = src_offset;
  wr->wc.dest_offset = dest_offset;
  wr->last_op = PAMI_OP_PUT_INITIATOR;
  wr->wc.op_complete = 1;
  wr->wc.local_send_complete = 2;
  pami_put_simple_t param;
  param.rma.dest                   = target_ep;
  param.rma.bytes                  = sizeof (nnti_bgqpami_work_completion);
  param.rma.cookie                 = (void *)&wr->wc.local_send_complete;
  param.rma.done_fn                = cb_done;
  param.addr.local      = (void *)&wr->wc;
  param.addr.remote     = (void *)src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.bgqpami.wc_addr;
  param.put.rdone_fn       = cb_done;


                    result = PAMI_Put(contexts[4], &param);

                       while (wr->wc.local_send_complete)
                          {
                                        result = PAMI_Context_trylock_advancev(&(contexts[4]), 1, 1000);
                          }

  nthread_unlock (&nnti_bgqpami_lock);
  nthread_lock (&target_buf_mem_hdl->wr_queue_lock);
  target_buf_mem_hdl->wr_queue.push_back (wr);
  nthread_unlock (&target_buf_mem_hdl->wr_queue_lock);
  insert_wr_wrhash (wr);
  log_debug (nnti_debug_level, "exit");

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
NNTI_bgqpami_wait (const NNTI_buffer_t * reg_buf,
		   const NNTI_buf_ops_t remote_op,
		   const int timeout, NNTI_status_t * status)
{
  NNTI_result_t
    nnti_rc = NNTI_OK;
  bgqpami_memory_handle *
    bgqpami_mem_hdl = NULL;
  pami_work_request *
    wr = NULL;
  bgqpami_connection *
    conn = NULL;
  int
    timeout_per_call;
  uint8_t
    retry_count = 0;


  log_debug (nnti_debug_level, "enter");

  bgqpami_mem_hdl = (bgqpami_memory_handle *) reg_buf->transport_private;

  if (timeout < 0)
    timeout_per_call = MIN_TIMEOUT;
  else
    timeout_per_call = (timeout < MIN_TIMEOUT) ? MIN_TIMEOUT : timeout;

  retry_count = 0;
  wr = bgqpami_mem_hdl->wr_queue.front ();
  bgqpami_mem_hdl->wr_queue.pop_front ();
  if ((wr->last_op == PAMI_OP_REGISTER_RDMA))
    {
      if ((remote_op == NNTI_SEND_SRC) || (remote_op == NNTI_GET_DST)
	  || (remote_op == NNTI_PUT_SRC))
	{
	  wr = bgqpami_mem_hdl->wr_queue.front ();
	  bgqpami_mem_hdl->wr_queue.pop_front ();
	}
    }
  if(wr == NULL){
	fprintf(stderr, "Process events no work request to work on\n");
	return NNTI_OK;
  }
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
	  //      del_wr_wrhash(wr);
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

      if ((bgqpami_mem_hdl->type == RDMA_TARGET_BUFFER)
	  || (bgqpami_mem_hdl->type == RECEIVE_BUFFER)
	  || (bgqpami_mem_hdl->type == GET_SRC_BUFFER)
	  || (bgqpami_mem_hdl->type == REQUEST_BUFFER)
	  || (bgqpami_mem_hdl->type == PUT_DST_BUFFER))
	{
	  bgqpami_mem_hdl->wr_queue.push_back (wr);
	}
      else
	{
	  del_wr_wrhash (wr);
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
NNTI_bgqpami_waitany (const NNTI_buffer_t ** buf_list,
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
	      assert (((bgqpami_memory_handle *) buf_list[i]->
		       transport_private)->type != REQUEST_BUFFER);
	    }
	}
    }
  assert (status);

  if (buf_count == 1)
    {
      nnti_rc = NNTI_bgqpami_wait (buf_list[0], remote_op, timeout, status);
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
		    NNTI_bgqpami_wait (buf_list[i], remote_op,
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

  log_debug (debug_level, "exit");

  trios_stop_timer ("NNTI_bgqpami_waitany", total_time);

  return (nnti_rc);
}

NNTI_result_t
NNTI_bgqpami_waitall (const NNTI_buffer_t ** buf_list,
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
NNTI_bgqpami_fini (const NNTI_transport_t * trans_hdl)
{
  /* Free up injection counter and global inj fifo and reception FIFO */
  return (NNTI_OK);
}

static NNTI_result_t
register_wc (pami_work_request * wr)
{
  pami_result_t 
    pami_result;
  size_t
    bytes_out;


/* This is for Work completion memory handle where ACK is received after PUT/GET */
  pami_result =
    PAMI_Memregion_create (contexts[3], &wr->wc,  sizeof (nnti_bgqpami_work_completion), &bytes_out, (pami_memregion_t *)wr->wc_mem_hdl);
  if (pami_result != PAMI_SUCCESS)
    {
      fprintf (stderr,
	       "PAMI memregion create failed in register memory Work completion handle\n");
      return (NNTI_EIO);
    }

  log_debug (nnti_debug_level, "exit  wr(%p)", wr);

  return ((NNTI_result_t) PAMI_SUCCESS);

}


static
  NNTI_result_t
register_memory (bgqpami_memory_handle * hdl, void *buf, uint64_t len,
		 const NNTI_buf_ops_t remote_op)
{
  NNTI_result_t
    rc = NNTI_OK;		/* return code */
  pami_result_t 
    pami_result;
  size_t
    bytes_out;

  assert (hdl);
  pami_result =
    PAMI_Memregion_create (contexts[3], buf, len, &bytes_out,  &hdl->mem_hdl);
  if (pami_result != PAMI_SUCCESS)
    fprintf (stderr, "PAMI memregion create failed in register memory \n");

  return (rc);
}

static int
unregister_memory (bgqpami_memory_handle * hdl)
{
  int
    rc = NNTI_OK;		/* return code */

  log_debug (nnti_debug_level, "exit  hdl(%p)", hdl);

  return (rc);
}


static
  NNTI_result_t
process_event (const NNTI_buffer_t * reg_buf,
	       const NNTI_buf_ops_t remote_op,
	       pami_work_request * wr, const int timeout)
{
  NNTI_result_t
    nnti_rc = NNTI_OK;
  bgqpami_memory_handle *
    bgqpami_mem_hdl = NULL;

  log_level
    debug_level = nnti_debug_level;

  bgqpami_mem_hdl = (bgqpami_memory_handle *) reg_buf->transport_private;

  log_debug (nnti_debug_level, "enter");
  long
    entry_time = trios_get_time_ms ();
  long
    elapsed_time = 0;

  debug_level = nnti_debug_level;
  switch (bgqpami_mem_hdl->type)
    {
    case SEND_BUFFER:
      while (wr->wc.remote_put_complete > 0)
	{
		pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[0]), 5, 1000);
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
      break;
    case PUT_SRC_BUFFER:
      wr->last_op = PAMI_OP_PUT_INITIATOR;
      if (wr->op_state == RDMA_WRITE_INIT)
	{
	  while (wr->wc.remote_put_complete > 0)
	    {
		pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[0]), 5, 1000);
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
	}
      break;
    case GET_DST_BUFFER:
      wr->last_op = PAMI_OP_GET_INITIATOR;
      if (wr->op_state == RDMA_READ_INIT)
	{
	  while (wr->wc.local_get_complete > 0)
	    {
		pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[1]), 5, 1000);
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
	}
      break;
    case REQUEST_BUFFER:
      {
	uint64_t
	  index = 0;
	bgqpami_request_queue_handle *
	  q = &transport_global_data.req_queue;

	wr->last_op = PAMI_OP_NEW_REQUEST;

	log_debug (debug_level,
		   "recv completion - reg_buf=%p current_req_index =%llu processing=%llu",
		   reg_buf, q->req_index, q->req_processed);

	while (q->req_index == q->req_processed)
	  {
/*
	    elapsed_time = (trios_get_time_ms () - entry_time);
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
		pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[2]), 1, 1000);

	  }
	index = q->req_processed;
	/* Shyamali wait for work completion buffer to showup */
	while (q->wc_buffer[q->req_processed].req_received != 1)
	  {
		pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[3]), 1, 1000);
	  }

	q->wc_buffer[q->req_processed].byte_offset = index * q->req_size;
	q->wc_buffer[q->req_processed].byte_len = q->req_size;
	wr->wc = q->wc_buffer[q->req_processed];
	wr->peer_rank = q->wc_buffer[q->req_processed].req_rank;
//	fprintf (stderr, "This request came from %d  process index %llu\n",
//		 wr->peer_rank, q->req_processed);
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
      wr->last_op = PAMI_OP_RECEIVE;
      log_debug (debug_level,
		 "receive buffer - recv completion - reg_buf==%p", reg_buf);

      while (wr->wc.op_complete != 1)
	{
		pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[4]), 1, 1000);
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
      wr->last_op = PAMI_OP_PUT_TARGET;
      while (wr->wc.op_complete != 1)
	{
		pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[0]), 5, 1000);

	}
      wr->peer_rank = wr->wc.req_rank;
      wr->op_state = RDMA_COMPLETE;
      break;
    case PUT_DST_BUFFER:
      wr->last_op = PAMI_OP_PUT_TARGET;
      while (wr->wc.op_complete != 1)
	{
	pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[2]), 2, 1000);
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
      break;
    case GET_SRC_BUFFER:
      wr->last_op = PAMI_OP_GET_TARGET;
      if (wr->op_state == RDMA_READ_INIT)
	{
	  while (wr->wc.op_complete != 1)
	    {
		pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[4]), 2, 1000);
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
      break;
    case RDMA_TARGET_BUFFER:
      if ((wr->last_op == PAMI_OP_GET_INITIATOR) ||
	  (wr->last_op == PAMI_OP_PUT_INITIATOR))
	{

	  if (wr->op_state == RDMA_TARGET_INIT)
	    {
	      while (wr->wc.op_complete != 1)
		{
		pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[2]), 3, 1000);
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
		pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[2]), 3, 1000);
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
  log_debug (nnti_debug_level, "exit");
  return (nnti_rc);
}


static void
create_peer (NNTI_peer_t * peer, int rank)
{
  sprintf (peer->url, "pami://%d/", rank);

  peer->peer.transport_id = NNTI_TRANSPORT_PAMI;
  peer->peer.NNTI_remote_process_t_u.bgqpami.pset_rank = rank;
}



/*
 * Loop over polling recv fifo 0  until everything arrives.
 */
static
  NNTI_result_t
recv_full (int rank, void *buf, size_t num)
{
  _recv_active = 1;
  while (_recv_active)
   pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[0]), 1, 1000);

  TRACE_ERR ((stderr,
	      "recv_full()  After advance  size recved %d  from (%zd)\n", num,
	      rank));
  memcpy (buf, &_recv_buffer, num);
  /* buf = &_recv_buffer; */
  _recv_active = 1;
  return NNTI_OK;
}


static
  NNTI_result_t
new_server_connection (bgqpami_connection * c, int rank)
{
  NNTI_result_t
    rc;
  bgqpami_request_queue_handle *
    q_hdl = &transport_global_data.req_queue;

  /*
   * Values passed through DMA MEMFIFO interface to permit initial connection.
   */
  struct
  {
    nnti_bgqpami_server_queue_attrs
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
	  sizeof (pami_memregion_t));

  sa_out.server_attrs.req_buffer_addr =
    (uint64_t) transport_global_data.req_queue.req_buffer;
  sa_out.server_attrs.req_size = transport_global_data.req_queue.req_size;
  sa_out.server_attrs.req_count = transport_global_data.req_queue.req_count;
  memcpy (&sa_out.server_attrs.req_mem_hdl, q_hdl->req_mem_hdl,
	  sizeof (pami_memregion_t));
  sa_out.server_attrs.wc_buffer_addr =
    (uint64_t) transport_global_data.req_queue.wc_buffer;
  memcpy (&sa_out.server_attrs.wc_mem_hdl,
	  &transport_global_data.req_queue.wc_mem_hdl,
	  sizeof (pami_memregion_t));
  sa_out.server_attrs.unblock_buffer_addr =
    (uint64_t) transport_global_data.req_queue.unblock_buffer;
  memcpy (&sa_out.server_attrs.unblock_mem_hdl,
	  &transport_global_data.req_queue.unblock_mem_hdl,
	  sizeof (pami_memregion_t));
   send_once ((char *)&sa_out, sizeof (sa_out), rank);

out:
  return rc;
}

static
  NNTI_result_t
insert_conn_rank (const uint32_t pset_rank, bgqpami_connection * conn)
{
  NNTI_result_t
    rc = NNTI_OK;

  nthread_lock (&nnti_conn_peer_lock);
  if (connections_by_rank.find (pset_rank) != connections_by_rank.end ())
    {
	conn = connections_by_rank[pset_rank];
    }
  connections_by_rank[pset_rank] = conn;	// add to connection map
  nthread_unlock (&nnti_conn_peer_lock);

  log_debug (nnti_debug_level, "peer connection added (conn=%p)", conn);

  return (rc);
}

static bgqpami_connection *
get_conn_rank (const uint32_t pset_rank)
{

  int
    key = pset_rank;
  bgqpami_connection *
    conn = NULL;

  nthread_lock (&nnti_conn_peer_lock);
  if (connections_by_rank.find (key) != connections_by_rank.end ())
    {
      conn = connections_by_rank[key];
    }

  nthread_unlock (&nnti_conn_peer_lock);

  if (conn != NULL)
    {
      log_debug (nnti_debug_level, "connection found");
      return conn;
    }

  log_debug (nnti_debug_level, "connection NOT found");
  return (NULL);
}

static bgqpami_connection *
del_conn_rank (const uint32_t pset_rank)
{
  bgqpami_connection *
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
      log_debug (nnti_debug_level, "connection found");
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

/*
 * At an explicit BYE message, or at finalize time, shut down a connection.
 * If descriptors are posted, defer and clean up the connection structures
 * later.
 */
static void
close_connection (bgqpami_connection * c)
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

  bgqpami_connection *
    conn = NULL;
  uint32_t
    peer_rank;
  int
    i;

  _recv_active = 1;
  while (_recv_active > 0)
    {
	   pami_result_t  result = PAMI_Context_trylock_advancev(&(contexts[1]), 1, 1000);
            if (result == PAMI_SUCCESS)
		  break;
    }
   	nthread_lock (&nnti_bgqpami_lock);
     	for(i=0; i<cur_conn_idx; i++) {
          	if (conn_list[i][1] == 0) { 
        		conn = get_conn_rank(conn_list[i][0]);
        		if (conn == NULL) {
				fprintf(stderr, "Connection not found check_poll for %d\n", i);
				break;
			}
			new_server_connection(conn, conn_list[i][0]);
			conn_list[i][1] = 1;
		}
        }
	cur_conn_idx = 0;
        nthread_unlock (&nnti_bgqpami_lock);
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
		     "error returned from nssi_bgqpami_server_listen_for_client: %d",
		     rc);
	  continue;
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
reset_req_index (bgqpami_request_queue_handle * req_queue_attrs)
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
fetch_server_req_buffer_offset (nnti_bgqpami_client_queue_attrs *
				local_req_queue_attrs,
				nnti_bgqpami_server_queue_attrs *
				remote_req_queue_attrs, uint64_t addend,
				uint64_t * prev_offset, int srv_rank)
{



  local_req_queue_attrs->req_index = -1;
  pami_endpoint_t target_ep;
  pami_result_t result = PAMI_Endpoint_create(client, (pami_task_t) srv_rank, 1, &target_ep);

  int active = 1;
  uint64_t  value;
  value  = 1;
  pami_rmw_t parameters;
  parameters.dest      = target_ep;
  //parameters.hints    = ;
    parameters.cookie    = &active;
    parameters.done_fn   = cb_done;
    parameters.local     = (void *)&local_req_queue_attrs->req_index;
    parameters.remote    = (void *)remote_req_queue_attrs->req_index_addr;
    parameters.value     = (void *)&value;
    parameters.operation = PAMI_ATOMIC_FETCH_ADD;
    parameters.type      = PAMI_TYPE_UNSIGNED_LONG_LONG;
   /* PAMI_ATOMIC_FETCH_ADD : local=remote and remote+=value */
                     uint64_t t0 = GetTimeBase();
 
   result = PAMI_Rmw(contexts[2], &parameters);
  
                           while (active)
                             {
                                          result = PAMI_Context_trylock_advancev(&(contexts[2]), 1, 1000);
                             }
  
                     uint64_t t1 = GetTimeBase();
                     uint64_t dt = t1-t0;
  
//   fprintf(stderr, " PAMI_Rmw local = %llu  in %llu cycles = %lf microseconds \n", local_req_queue_attrs->req_index,  (long long unsigned) dt, dt/1600.0 ); 

  *prev_offset = local_req_queue_attrs->req_index;
  return (0);
}

static int
send_req (nnti_bgqpami_client_queue_attrs * local_req_queue_attrs,
	  nnti_bgqpami_server_queue_attrs * remote_req_queue_attrs,
	  uint64_t offset, const NNTI_buffer_t * reg_buf, int rank,
	  pami_work_request * wr)
{
  bgqpami_memory_handle *
    src_buf_mem_hdl = NULL;
  src_buf_mem_hdl = (bgqpami_memory_handle *) reg_buf->transport_private;
  uint64_t
    length;
  pami_endpoint_t target_ep;
  length = reg_buf->payload_size;
    wr->wc.remote_put_complete = 2;	/* callback will decr  when recv done on remote side */
    wr->wc.local_send_complete = 1;	/* First callback returns after  injected into local TORUS */
    pami_result_t result = PAMI_Endpoint_create(client, (pami_task_t) rank, 1, &target_ep);
  pami_put_simple_t parameters;
  parameters.rma.dest                   = target_ep;
  parameters.rma.bytes                  = length;
  parameters.rma.cookie                 = (void *)&wr->wc.remote_put_complete;
  parameters.rma.done_fn                = decrement;
  parameters.addr.local      = (void *)reg_buf->payload;
  parameters.addr.remote     = (void *)remote_req_queue_attrs->req_buffer_addr + offset; 
  parameters.put.rdone_fn       = decrement;

                   uint64_t t0 = GetTimeBase();
//  fprintf(stderr, "buffer address %p remote address %p offset %llu\n", parameters.addr.local, parameters.addr.remote, offset);

                    result = PAMI_Put(contexts[3], &parameters);
                       while (wr->wc.remote_put_complete)
                          {
                                        result = PAMI_Context_trylock_advancev(&(contexts[3]), 1, 1000);
                          }

                                               uint64_t t1 = GetTimeBase();
                                                 uint64_t dt = t1-t0;
  //fprintf(stderr, " send req  done %llu cycles = %lf microseconds\n", (long long unsigned)dt, dt/1600.0);

  return (0);
}

static int
send_req_wc (nnti_bgqpami_client_queue_attrs * client_q,
	     nnti_bgqpami_server_queue_attrs * server_q,
	     const NNTI_buffer_t * reg_buf, uint64_t offset,
	     pami_work_request * wr, int rank)
{
  bgqpami_memory_handle *
    bgqpami_mem_hdl = NULL;

  bgqpami_mem_hdl = (bgqpami_memory_handle *) reg_buf->transport_private;
  pami_endpoint_t target_ep;
    wr->wc.remote_put_complete = 2;     /* callback will decr  when recv done on remote side */
    wr->wc.local_send_complete = 1;     /* First callback returns after  injected into local TORUS */
    pami_result_t result = PAMI_Endpoint_create(client, (pami_task_t) rank, 1, &target_ep);
  pami_put_simple_t parameters;
  parameters.rma.dest                   = target_ep;
  parameters.rma.bytes                  = sizeof (nnti_bgqpami_work_completion);
  parameters.rma.cookie                 = (void *)&wr->wc.remote_put_complete;
  parameters.rma.done_fn                = cb_done;
  parameters.addr.local      = (void *)&wr->wc;
  parameters.addr.remote     = (void *)server_q->wc_buffer_addr + offset;
  parameters.put.rdone_fn       = cb_done;

                   uint64_t t0 = GetTimeBase();

                    result = PAMI_Put(contexts[3], &parameters);

                       while (wr->wc.remote_put_complete)
                          {
                                        result = PAMI_Context_trylock_advancev(&(contexts[3]), 1, 1000);
                          }

//	fprintf(stderr, "Send request wc completed \n"); 
  return (0);

}



static int
request_send (nnti_bgqpami_client_queue_attrs * client_q,
	      nnti_bgqpami_server_queue_attrs * server_q,
	      const NNTI_buffer_t * reg_buf, int rank)
{
  int
    rc = 0;
  bgqpami_memory_handle *
    bgqpami_mem_hdl = NULL;
  uint64_t
    offset;
  pami_work_request *
    wr = NULL;

  uint32_t
    wc_size = sizeof (nnti_bgqpami_work_completion);

  bgqpami_mem_hdl = (bgqpami_memory_handle *) reg_buf->transport_private;

  assert (bgqpami_mem_hdl);
  wr = (pami_work_request *) malloc (sizeof (pami_work_request));
  memset (wr, 0, sizeof (pami_work_request));
  assert (wr);

  wr->reg_buf = reg_buf;

  bgqpami_mem_hdl->wr_queue.push_back (wr);

  insert_wr_wrhash (wr);

  log_debug (nnti_debug_level, "enter");

  log_debug (nnti_debug_level, "calling fetch_add_buffer_offset()");
  rc = fetch_server_req_buffer_offset (client_q, server_q, 1, &offset, rank);

  log_debug (nnti_debug_level, "calling send_req()");
  rc =
    send_req (client_q, server_q, offset * server_q->req_size, reg_buf, rank,
	      wr);
  if (rc != NNTI_OK)
    log_error (nnti_debug_level, "send_req() failed: %d", rc);
  wr->wc.op = PAMI_OP_SEND_REQ;
  wr->wc.req_received = 1;
  wr->wc.req_rank = transport_global_data.myrank;
  wr->wc.byte_len = reg_buf->payload_size;
  wr->wc.byte_offset = client_q->req_index * server_q->req_size;
  wr->wc.src_offset = 0;
  wr->wc.dest_offset = client_q->req_index * server_q->req_size;
  wr->peer_rank = rank;

/* Shyamali */ 
  log_debug (nnti_debug_level, "calling send_wc()");
  rc = send_req_wc (client_q, server_q, reg_buf, offset * wc_size, wr, rank);
 /* Shyamali END */
  log_debug (nnti_debug_level, "exit");
  return (0);
}



static int
client_req_queue_init (bgqpami_connection * c)
{
  int
    rc = 0;
  size_t
    bytes;

  nnti_bgqpami_client_queue_attrs *
    q = &c->queue_local_attrs;
  q->req_index = 0;
  q->req_index_addr = (uint64_t) & q->req_index;
  PAMI_Memregion_create (&q->req_index_mem_hdl, 
			  (void *) q->req_index_addr, sizeof (uint64_t), &bytes,
                          (pami_memregion_t *)&q->req_index_mem_hdl);
  return (rc);
}

static int
client_req_queue_destroy (bgqpami_connection * c)
{
  return (0);
}

static int
server_req_queue_init (bgqpami_request_queue_handle * q,
		       char *buffer, uint64_t req_size, uint64_t req_count)
{
  int
    i;
  size_t
    bytes;
  uint64_t
    unblock_size;
  bgqpami_memory_handle *
    bgqpami_mem_hdl = (bgqpami_memory_handle *) q->reg_buf->transport_private;


  q->req_buffer = buffer;
  q->req_size = req_size;
  q->req_count = req_count;
  q->req_buffer_size = req_count * req_size;

  q->wc_buffer_size = q->req_count * sizeof (nnti_bgqpami_work_completion);
  q->wc_buffer =
    (nnti_bgqpami_work_completion *) calloc (q->req_count,
					     sizeof
					     (nnti_bgqpami_work_completion));

  unblock_size = MAX_CONNECTION * sizeof (nnti_bgqpami_work_completion);
  q->unblock_buffer =
    (nnti_bgqpami_fetch_index_req *) calloc (MAX_CONNECTION,
					     sizeof
					     (nnti_bgqpami_fetch_index_req));
  PAMI_Memregion_create (contexts[2], buffer,  q->req_buffer_size, &bytes, &q->req_mem_hdl);
  PAMI_Memregion_create (contexts[2], q->wc_buffer,  q->wc_buffer_size,&bytes, &q->wc_mem_hdl);
  PAMI_Memregion_create (contexts[2],q->unblock_buffer, unblock_size, &bytes, &q->unblock_mem_hdl);
  q->req_index = 0;
  q->req_index_addr = (uint64_t) & q->req_index;

  q->req_processed = 0;
  q->req_processed_reset_limit = q->req_count;
  bgqpami_mem_hdl->type = REQUEST_BUFFER;
  PAMI_Memregion_create (contexts[2], (void *) q->req_index_addr,  sizeof(uint64_t), &bytes, &q->req_index_mem_hdl);
  for (i = 0; i < MAX_CONNECTION; ++i)
    q->unblock_buffer[i].req_rank = -1;

  start_connection_listener_thread ();
  return 0;
}

static int
server_req_queue_destroy (bgqpami_request_queue_handle * q)
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
    hash6432shift ((uint64_t) buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.
		   buf);

  nthread_lock (&nnti_buf_bufhash_lock);
  assert (buffers_by_bufhash.find (h) == buffers_by_bufhash.end ());
  buffers_by_bufhash[h] = buf;
  nthread_unlock (&nnti_buf_bufhash_lock);

  log_debug (nnti_debug_level, "bufhash buffer added (buf=%p)", buf);

  return (rc);
}

static NNTI_buffer_t *
get_buf_bufhash (const uint32_t bufhash)
{
  NNTI_buffer_t *
    buf = NULL;

  log_debug (nnti_debug_level, "looking for bufhash=%llu",
	     (uint64_t) bufhash);
  nthread_lock (&nnti_buf_bufhash_lock);
  if (buffers_by_bufhash.find (bufhash) != buffers_by_bufhash.end ())
    {
      buf = buffers_by_bufhash[bufhash];
    }
  nthread_unlock (&nnti_buf_bufhash_lock);

  if (buf != NULL)
    {
      log_debug (nnti_debug_level, "buffer found (buf=%p)", buf);
      return buf;
    }

  log_debug (nnti_debug_level, "buffer NOT found");
  print_bufhash_map ();

  return (NULL);
}

static NNTI_buffer_t *
del_buf_bufhash (NNTI_buffer_t * buf)
{
  uint32_t
    h =
    hash6432shift ((uint64_t) buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.
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
insert_wr_wrhash (pami_work_request * wr)
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

static pami_work_request *
get_wr_wrhash (const uint32_t wrhash)
{
  pami_work_request *
    wr = NULL;

  log_debug (nnti_debug_level, "looking for wrhash=%llu", (uint64_t) wrhash);
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

static pami_work_request *
del_wr_wrhash (pami_work_request * wr)
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
      log_debug (debug_level, "work request found");
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
  bgqpami_memory_handle *
    pami_mem_hdl = NULL;

  log_debug (nnti_debug_level, "enter");

  pami_mem_hdl = (bgqpami_memory_handle *) reg_buf->transport_private;
  assert (pami_mem_hdl);

  if (pami_mem_hdl->wr_queue.empty ())
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
  bgqpami_memory_handle *
    pami_mem_hdl = NULL;
  pami_work_request *
    wr = NULL;

  log_debug (nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

  pami_mem_hdl = (bgqpami_memory_handle *) reg_buf->transport_private;

  if (is_wr_queue_empty (reg_buf) == TRUE)
    {
      log_debug (nnti_debug_level,
		 "work request queue is empty - return FALSE");
      rc = FALSE;
    }
  else
    {
      wr = pami_mem_hdl->wr_queue.front ();
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
is_wr_complete (pami_work_request * wr)
{
  int8_t
    rc = FALSE;
  bgqpami_memory_handle *
    pami_mem_hdl = NULL;

  pami_mem_hdl = (bgqpami_memory_handle *) wr->reg_buf->transport_private;
  assert (pami_mem_hdl);

  switch (pami_mem_hdl->type)
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
repost_recv_work_request (NNTI_buffer_t * reg_buf, pami_work_request * wr)
{
  bgqpami_memory_handle *
    bgqpami_mem_hdl = NULL;

  log_debug (nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

  bgqpami_mem_hdl = (bgqpami_memory_handle *) reg_buf->transport_private;
  assert (bgqpami_mem_hdl);

  memset (&wr->wc, 0, sizeof (nnti_bgqpami_work_completion));

  wr->last_op = 0;

  memset (&wr->op_state, 0, sizeof (bgqpami_op_state_t));
  register_wc (wr);

  bgqpami_mem_hdl->wr_queue.push_back (wr);

  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.wc_addr =
    (uint64_t) & wr->wc;
  reg_buf->buffer_addr.NNTI_remote_addr_t_u.bgqpami.wc_mem_hdl = (uint64_t)wr->wc_mem_hdl;
  log_debug (nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

  return (NNTI_OK);
}
