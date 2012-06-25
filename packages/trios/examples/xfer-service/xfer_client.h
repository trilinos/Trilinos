/*
 * xfer_client.h
 *
 *  Created on: Aug 22, 2011
 *      Author: raoldfi
 */

#ifndef XFER_CLIENT_H_
#define XFER_CLIENT_H_

#include <string>
#include <limits.h>

#include "Trios_logger.h"



enum IO_METHODS {
    XFER_WRITE_ENCODE_SYNC = 0,
    XFER_WRITE_ENCODE_ASYNC,
    XFER_WRITE_RDMA_SYNC,
    XFER_WRITE_RDMA_ASYNC,
    XFER_READ_ENCODE_SYNC,
    XFER_READ_ENCODE_ASYNC,
    XFER_READ_RDMA_SYNC,
    XFER_READ_RDMA_ASYNC
};

enum DIST_SCHEME {
    XFER_BLOCK_DISTRIBUTION,
    XFER_ROUND_ROBIN_DISTRIBUTION
};

enum MPI_IO_METHODS {
    XFER_MPI_SEND=0,
    XFER_MPI_ISEND,
    XFER_MPI_RECV,
    XFER_MPI_IRECV,
    XFER_MPI_PUT,
    XFER_MPI_GET
};


/**
 * Options and arguments passed to the client driver.
 */
struct xfer_args {
        bool client_flag;
        bool server_flag;
        int num_servers;  // used for the server only
        int num_threads;  // used for the server only
        bool block_distribution; // how to assign clients to servers
        std::string server_url;
        int transport;
        std::string transport_name;
        int len;
        int io_method;
        std::string url_file;
        std::string io_method_name;
        log_level debug_level;
        std::string logfile;
        int num_trials;
        int num_reqs;
        std::string result_file;
        std::string result_file_mode;
        int timeout;
        int delay;
        int num_retries;
        bool validate_flag;
};


#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)



#else /* K&R C */
#endif




#ifdef __cplusplus
}
#endif


#endif /* XFER_CLIENT_H_ */
