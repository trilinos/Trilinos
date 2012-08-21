/*-------------------------------------------------------------------------*/
/**  @file nnti_xdr.x
 *
 *   @brief XDR definitions for types used in the NNTI transport interface.
 *
 *   @author Todd Kordenbrock (thkorde\@sandia.gov).
 *
 */

#include "Trios_config.h"

#ifdef RPC_HDR
%#include <rpc/types.h>
%#include <rpc/xdr.h>
%#include "Trios_config.h"
%#include <sys/param.h>
%#include <errno.h>
%#include <stdint.h>

#endif

#ifdef RPC_XDR
%#include "Trios_nnti_xdr.h"
%#include "Trios_config.h"
%#ifdef HAVE_TRIOS_XDR_U_INT16_T
%#define xdr_uint16_t xdr_u_int16_t
%#endif
%#ifdef HAVE_TRIOS_XDR_U_INT32_T
%#define xdr_uint32_t xdr_u_int32_t
%#endif
%#ifdef HAVE_TRIOS_XDR_U_INT64_T
%#define xdr_uint64_t xdr_u_int64_t
%#endif

#endif


/** @addtogroup nnti_types
 *  @{
 */

#ifdef HAVE_TRIOS_U_INT32_T
typedef u_int32_t uint32_t
#endif

/**
 * @brief The <tt>\ref nssi_size</tt> type is used to for ``size'' variables.
 */
typedef uint64_t NNTI_size;
typedef uint64_t NNTI_ssize;

/**
 * @brief Length of a hostname
 */
const NNTI_HOSTNAME_LEN = 32;

/**
 * @brief Length of a URL
 */
const NNTI_URL_LEN = 128;

/**
 * @brief The size of a buffer used for receiving requests (should be configurable via URL).
 */
const NNTI_REQUEST_BUFFER_SIZE=1200;
/**
 * @brief The size of a buffer used for receiving results (should be configurable via URL).
 */
const NNTI_RESULT_BUFFER_SIZE=1200;


/**
 * @brief Enumerator of the transport mechanisms supported by NNTI.
 *
 * The <tt>\ref NNTI_transport_id_t</tt> enumerator provides integer values
 * to represent the supported transport mechanisms.
 */
enum NNTI_transport_id_t {
    /** @brief Use Portals to transfer rpc requests. */
    NNTI_TRANSPORT_PORTALS,

    /** @brief Use Infiniband to transfer rpc requests. */
    NNTI_TRANSPORT_IB,

    /** @brief Use Cray Gemini to transfer rpc requests. */
    NNTI_TRANSPORT_GEMINI,

    /** @brief Use Cray Gemini to transfer rpc requests. */
    NNTI_TRANSPORT_MPI,

    /** @brief Use Cray LUC to transfer rpc requests. */
    NNTI_TRANSPORT_LUC,

    /** @brief Use a local buffer (no remote operations). */
    NNTI_TRANSPORT_LOCAL,

    /** @brief No operations permitted. */
    NNTI_TRANSPORT_NULL
};

/**
 * @brief The number of transport mechanisms supported by NNTI.
 */
const NNTI_TRANSPORT_COUNT = 7;


/**
 * @brief Enumerator of results that NNTI functions could generate.
 *
 * The <tt>\ref NNTI_result_t</tt> enumerator provides integer values
 * for NNTI function outcomes.
 */
enum NNTI_result_t {
    /** @brief The function completed successfully. */
    NNTI_OK = 0,

    /** @brief Unspecified I/O error. */
    NNTI_EIO = 1001,

    /** @brief The size of the message is larger than the supported maximum size. */
    NNTI_EMSGSIZE,

    /** @brief The operation or process has been canceled. */
    NNTI_ECANCELED,

    /** @brief An operation timed out. */
    NNTI_ETIMEDOUT,

    /** @brief The value of the variable or parameter is invalid. */
    NNTI_EINVAL,

    /** @brief  No memory is available. */
    NNTI_ENOMEM,

    /** @brief No such entry. */
    NNTI_ENOENT,

    /** @brief Unsupported operation. */
    NNTI_ENOTSUP,

    /** @brief The item already exists. */
    NNTI_EEXIST,

    /** @brief Unsuccessful RPC operation. */
    NNTI_EBADRPC,

    /** @brief Not initialized. */
    NNTI_ENOTINIT,

    /** @brief Insufficient priveleges to perform operation. */
    NNTI_EPERM,

    /** @brief An operation would have blocked. */
    NNTI_EWOULDBLOCK,

    /** @brief Operation was interupted, but possibly recoverable. */
    NNTI_EAGAIN

};



/**
 * @brief Network node identifier (used by Portals).
 *
 * The <tt>\ref NNTI_nid</tt> type identifies a particular node.
 */
typedef uint32_t NNTI_nid;

/**
 * @brief Process identifier (used by Portals).
 *
 * The <tt>\ref NNTI_pid</tt> type identifies a particular process.
 */
typedef uint32_t NNTI_pid;



/**
 * @brief Remote process identifier for Portals.
 *
 * The <tt>\ref NNTI_portals_process_t</tt> identifies a particular process
 * on a particular node.
 */
struct NNTI_portals_process_t {
    /** @brief Portals Node ID. */
    NNTI_nid nid;
    /** @brief Portals Process ID. */
    NNTI_pid pid;
};




/**
 * @brief Binary encoding of a TCP/IP host address.
 *
 * The <tt>\ref NNTI_ip_addr</tt> type identifies a particular node.
 */
typedef uint32_t NNTI_ip_addr;
/**
 * @brief TCP port in NBO.
 *
 * The <tt>\ref NNTI_tcp_addr</tt> type identifies a particular port.
 */
typedef uint16_t NNTI_tcp_port;
/**
 * @brief The queue pair number of an IB connection.
 *
 * The <tt>\ref NNTI_qp_num</tt> type identifies a particular connection
 * between two IB hosts.
 */
typedef uint32_t NNTI_qp_num;

/**
 * @brief Remote process identifier for IB.
 *
 * The <tt>\ref NNTI_ib_process_t</tt> identifies a particular process
 * on a particular node.  If a connection has been established to the
 * represented process, then that connection is identified by 'qp_num'.
 */
struct NNTI_ib_process_t {
/*    /** @brief name on the peer */
/*    string        hostname<NSSI_HOSTNAME_LEN>; */
    /** @brief IP address encoded in Network Byte Order */
    NNTI_ip_addr  addr;
    /** @brief TCP port encoded in Network Byte Order */
    NNTI_tcp_port port;
    /** @brief IB connection queue pair number */
    NNTI_qp_num   qpn;
};







/**
 * @brief The ID of a Cray LUC Endpoint.
 *
 * The <tt>\ref NNTI_luc_endpoint_id</tt> type identifies a particular LUC process
 * on a particular node.
 */
typedef uint64_t NNTI_luc_endpoint_id;    /* id used to contact a LUC program */

/**
 * @brief Remote process identifier for Cray LUC.
 *
 * The <tt>\ref NNTI_luc_process_t</tt> identifies a particular process
 * on a particular node.
 */
struct NNTI_luc_process_t {
    /** @brief ID used to contact a LUC program */
    NNTI_luc_endpoint_id req_endpoint_id;
    /** @brief ID used to do RDMA transfers with a LUC program */
    NNTI_luc_endpoint_id rdma_endpoint_id;
};










/**
 * @brief The instance ID of a Gemini process.
 *
 * The <tt>\ref NNTI_inst_id</tt> type identifies a particular process
 * within a communication domain.
 */
typedef uint32_t NNTI_instance_id;

/**
 * @brief Remote process identifier for Gemini.
 *
 * The <tt>\ref NNTI_gni_process_t</tt> identifies a particular process
 * on a particular node.  If a connection has been established to the
 * represented process, then that connection is identified by 'inst_id'.
 */
struct NNTI_gni_process_t {
    /** @brief IP address encoded in Network Byte Order */
    NNTI_ip_addr  addr;
    /** @brief TCP port encoded in Network Byte Order */
    NNTI_tcp_port port;
    /** @brief Gemini process instance ID */
    NNTI_instance_id  inst_id;
};






/**
 * @brief Remote process identifier for MPI.
 *
 * The <tt>\ref NNTI_mpi_process_t</tt> identifies a particular process
 * on a particular node.
 */
struct NNTI_mpi_process_t {
    /** @brief MPI rank. */
    int rank;
};





/**
 * @brief A structure to represent a remote processes.
 *
 * The <tt>NNTI_remote_process_t</tt> structure contains the
 * transport specific info needed to identify a process running
 * on a remote node.
 */
#if defined(RPC_HDR) || defined(RPC_XDR)
union NNTI_remote_process_t switch (NNTI_transport_id_t transport_id) {
    /** @brief The Portals representation of a process on the network. */
    case NNTI_TRANSPORT_PORTALS: NNTI_portals_process_t portals;
    /** @brief The IB representation of a process on the network. */
    case NNTI_TRANSPORT_IB:      NNTI_ib_process_t      ib;
    /** @brief The Cray LUC representation of a process on the network. */
    case NNTI_TRANSPORT_LUC:     NNTI_luc_process_t     luc;
    /** @brief The Cray Gemini representation of a process on the network. */
    case NNTI_TRANSPORT_GEMINI:  NNTI_gni_process_t     gni;
    /** @brief The MPI representation of a process on the network. */
    case NNTI_TRANSPORT_MPI:     NNTI_mpi_process_t     mpi;
};
#else
union NNTI_remote_process_t {
    /** @brief The Portals representation of a process on the network. */
    NNTI_portals_process_t portals;
    /** @brief The IB representation of a process on the network. */
    NNTI_ib_process_t      ib;
    /** @brief The Cray LUC representation of a process on the network. */
    NNTI_luc_process_t     luc;
    /** @brief The Cray Gemini representation of a process on the network. */
    NNTI_gni_process_t     gni;
    /** @brief The MPI representation of a process on the network. */
    NNTI_mpi_process_t     mpi;
};
#endif


/**
 * @brief Handle to an NNTI process.
 *
 * This is the datatype used by NNTI clients to reference another process.
 * Use this handle to move data to/from the process.
 */
struct NNTI_peer_t {
    /** @brief string encoding of a process on the network */
    opaque url[NNTI_URL_LEN];

    /** @brief binary encoding of a process on the network */
    NNTI_remote_process_t peer;
};






/**
 * @brief Definition for match bits in Portals.
 */
typedef uint64_t NNTI_match_bits;

/**
 * @brief Definition for a remote buffer id.
 */
enum NNTI_portals_indices {
    NNTI_REQ_PT_INDEX = 1,   /* where to send requests */
    NNTI_RECV_PT_INDEX,      /* where to receive data */
    NNTI_DATA_PT_INDEX,      /* where to put/get data */
    NNTI_LONG_ARGS_PT_INDEX, /* where to fetch long args */
    NNTI_LONG_RES_PT_INDEX   /* where to fetch long results */
};


/**
 * @brief RDMA address used for the Portals implementation.
 */
struct NNTI_portals_rdma_addr_t {
    /** @brief The memory buffer ID. */
    NNTI_portals_indices buffer_id;
    /** @brief The match bits (required by Portals). */
    NNTI_match_bits      match_bits;
    /** @brief Size of the the memory buffer. */
    uint32_t size;
};










/**
 * @brief RDMA address used for the InfiniBand implementation.
 */
struct NNTI_ib_rdma_addr_t {
    /** @brief Address of the memory buffer cast to a uint64_t. */
    uint64_t buf;
    /** @brief The key that a remote processes needs to access this buffer. */
    uint32_t key;
    /** @brief Size of the the memory buffer. */
    uint32_t size;

    /** @brief Address of the ACK buffer cast to a uint64_t. */
    uint64_t ack_buf;
    /** @brief The key that a remote processes needs to access the ACK buffer. */
    uint32_t ack_key;
    /** @brief Size of the the ACK buffer. */
    uint32_t ack_size;
};




/**
 * @brief RDMA address used for LUC implementation.
 */
struct NNTI_luc_rdma_addr_t {
    /** @brief Address of the memory buffer cast to a uint64_t. */
    uint64_t buf;
    /** @brief Size of the the memory buffer. */
    uint32_t size;
};

enum NNTI_gni_buffer_type_t {
    NNTI_GNI_RDMA_INITIATOR,
    NNTI_GNI_RDMA_TARGET,
    NNTI_GNI_SEND_SRC,
    NNTI_GNI_REQUEST_BUFFER
};

struct NNTI_gni_mem_hdl_t {
    uint64_t qword1;
    uint64_t qword2;
};

/**
 * @brief RDMA address used for the Gemini implementation.
 */
struct NNTI_gni_rdma_addr_t {
    /** @brief Address of the memory buffer cast to a uint64_t. */
    uint64_t buf;
    /** @brief Size of the the memory buffer. */
    uint32_t size;
    /** @brief The key that a remote processes needs to access this buffer. */
    NNTI_gni_mem_hdl_t mem_hdl;

    /** @brief Address of the Work Completion buffer cast to a uint64_t. */
    uint64_t wc_addr;
    /** @brief The key that a remote processes needs to access the Work Completion buffer. */
    NNTI_gni_mem_hdl_t wc_mem_hdl;

    /** @brief Identifies the type of buffer (FMA Short Msg or RDMA) */
    NNTI_gni_buffer_type_t type;
};



/**
 * @brief RDMA address used for the MPI implementation.
 */
struct NNTI_mpi_rdma_addr_t {
    /** @brief The MPI tag for RTR msg. */
    NNTI_match_bits rtr_tag;
    /** @brief The MPI tag for RTS msg. */
    NNTI_match_bits rts_tag;
    /** @brief The MPI tag for data msg. */
    NNTI_match_bits data_tag;
    /** @brief Size of the the memory buffer. */
    uint32_t        size;
};




/**
 * @brief A structure to represent a remote memory region.
 *
 * The <tt>NNTI_remote_addr_t</tt> structure contains the
 * transport specific info needed to identify a memory region
 * on a remote node.
 */
#if defined(RPC_HDR) || defined(RPC_XDR)
union NNTI_remote_addr_t switch (NNTI_transport_id_t transport_id) {
    /** @brief The Portals representation of a memory region. */
    case NNTI_TRANSPORT_PORTALS: NNTI_portals_rdma_addr_t portals;
    /** @brief The IB representation of a memory region. */
    case NNTI_TRANSPORT_IB:      NNTI_ib_rdma_addr_t      ib;
    /** @brief The Cray LUC representation of a memory region. */
    case NNTI_TRANSPORT_LUC:     NNTI_luc_rdma_addr_t     luc;
    /** @brief The Cray Gemini representation of a memory region. */
    case NNTI_TRANSPORT_GEMINI:  NNTI_gni_rdma_addr_t     gni;
    /** @brief The MPI representation of a memory region. */
    case NNTI_TRANSPORT_MPI:     NNTI_mpi_rdma_addr_t     mpi;
};
#else
union NNTI_remote_addr_t {
    /** @brief The Portals representation of a memory region. */
    NNTI_portals_rdma_addr_t portals;
    /** @brief The IB representation of a memory region. */
    NNTI_ib_rdma_addr_t      ib;
    /** @brief The Cray LUC representation of a memory region. */
    NNTI_luc_rdma_addr_t     luc;
    /** @brief The Cray Gemini representation of a memory region. */
    NNTI_gni_rdma_addr_t     gni;
    /** @brief The MPI representation of a memory region. */
    NNTI_mpi_rdma_addr_t     mpi;
};
#endif


/**
 * @brief The operations that are permitted on a buffer
 *
 * The <tt>NNTI_buf_ops_t</tt> enum defines the operations that can be
 * performed on a buffer.  These operations can be OR'd together to create
 * a multipurpose buffer.
 */
enum NNTI_buf_ops_t {
    /** @brief this buffer can be put from */
    NNTI_PUT_SRC=1,
    /** @brief this buffer can be put into */
    NNTI_PUT_DST=2,
    /** @brief this buffer can be got from */
    NNTI_GET_SRC=4,
    /** @brief this buffer can be got into */
    NNTI_GET_DST=8,
    /** @brief this buffer can be sent from */
    NNTI_SEND_SRC=16,
    /** @brief this buffer can be received into */
    NNTI_RECV_DST=32,
    /** @brief this buffer has multiple receive slots */
    NNTI_RECV_QUEUE=64
};


/**
 * @brief handle to a memory buffer prepared by NNTI_register_memory
 *
 * The <tt>NNTI_buffer_t</tt> structure contains the
 * location of a buffer on the network.  This is all the info
 *  a peer needs to put/get this buffer.
 */
struct NNTI_buffer_t {
    /** @brief the transport where this buffer is registered */
    NNTI_transport_id_t transport_id;

    /** @brief the process in which this buffer resides */
    NNTI_peer_t        buffer_owner;
    /** @brief the remote address at which this buffer resides */
    NNTI_remote_addr_t buffer_addr;

    /** @brief permitted operations */
    NNTI_buf_ops_t ops;

    /** @brief Size of this buffer. */
    uint64_t     payload_size;
    /** @brief Local address of the memory buffer cast to a uint64_t. */
    uint64_t     payload;

    /** @brief Private storage (cast to a uint64_t). */
    uint64_t     transport_private;
};

/**
 * @brief A header that is prepended to all requests and results.
 *
 */
struct NNTI_transport_header_t {
    /** @brief a piece of magic to describe the type of message contained in the payload */
    uint32_t message_magic;
    /** @brief some identifier the transport can use to create a thread of messages */
    uint32_t exchange_id;
};


/**
 * @brief The status of an NNTI operation.
 *
 * The <tt>NNTI_status_t</tt> structure contains the result of an NNTI buffer
 * operation.  The result and op fields will always contain reasonable values.
 * If the operation was successful (result==NNTI_OK), then the remaining
 * fields will also contain reasonable values.  If the operation was
 * unsuccessful, then the remaining fields are undefined.
 */
struct NNTI_status_t {
    /** @brief The operation performed on the buffer. */
    NNTI_buf_ops_t op;
    /** @brief The result code for this operation. */
    NNTI_result_t  result;

    /** @brief The address of the local buffer used by this operation. */
    uint64_t start;
    /** @brief The offset into the local buffer used by this operation. */
    uint64_t offset;
    /** @brief The number of bytes used by this operation. */
    uint64_t length;

    /** @brief The peer that was the data source for this operation. */
    NNTI_peer_t src;
    /** @brief The peer that was the data destination for this operation. */
    NNTI_peer_t dest;
};


/**
 * @brief The external representation of a configured transport.
 */
struct NNTI_transport_t {
    /** @brief The transport id. */
    NNTI_transport_id_t  id;
    /** @brief A reference to my process that can be sent to a peer so the peer can contact me. */
    NNTI_peer_t          me;
};


/** @} */
