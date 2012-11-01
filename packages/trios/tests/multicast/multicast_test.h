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
 * multicast_client.h
 *
 *  Created on: Nov 14, 2011
 *      Author: thkorde
 */

#ifndef MULTICAST_TEST_H_
#define MULTICAST_TEST_H_

#include <string>
#include <limits.h>

#include "Trios_logger.h"



enum IO_METHODS {
    MULTICAST_EMPTY_REQUEST_SYNC,
    MULTICAST_EMPTY_REQUEST_ASYNC,
    MULTICAST_GET_SYNC,
    MULTICAST_GET_ASYNC,
    MULTICAST_PUT_SYNC,
    MULTICAST_PUT_ASYNC
};


/**
 * Options and arguments passed to the client driver.
 */
struct multicast_args {
        bool client_flag;
        bool server_flag;
        int transport;
        std::string transport_name;
        int len;
        int io_method;
        std::string server_url[2];
        std::string url_file[2];
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


#endif /* MULTICAST_TEST_H_ */
