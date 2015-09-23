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

/*
 * netperf_client.h
 *
 *  Created on: Aug 22, 2011
 *      Author: raoldfi
 */

#ifndef NETPERF_CLIENT_H_
#define NETPERF_CLIENT_H_

#include <string>
#include <limits.h>

#include "Trios_logger.h"



enum IO_METHODS {
    NETPERF_WRITE_RDMA_SYNC=0,
    NETPERF_WRITE_RDMA_ASYNC,
    NETPERF_READ_RDMA_SYNC,
    NETPERF_READ_RDMA_ASYNC
};

enum DIST_SCHEME {
    NETPERF_BLOCK_DISTRIBUTION,
    NETPERF_ROUND_ROBIN_DISTRIBUTION
};

enum MPI_IO_METHODS {
    NETPERF_MPI_SEND=0,
    NETPERF_MPI_ISEND,
    NETPERF_MPI_RECV,
    NETPERF_MPI_IRECV,
    NETPERF_MPI_PUT,
    NETPERF_MPI_GET
};


/**
 * Options and arguments passed to the client driver.
 */
struct netperf_args {
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


#endif /* NETPERF_CLIENT_H_ */
