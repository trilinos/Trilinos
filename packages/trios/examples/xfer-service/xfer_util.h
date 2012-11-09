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
 * xfer_util.h
 *
 *  Created on: May 31, 2012
 *      Author: raoldfi
 */

#ifndef XFER_UTIL_H_
#define XFER_UTIL_H_

#include <mpi.h>
#include <string>
#include <vector>
#include "xfer_service_args.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int xfer_write_server_url_file(
        std::string url_fname,
        std::string my_url,
        MPI_Comm comm);


extern void xfer_read_server_url_file(
        const char *path,
        std::vector<std::string> &urlbuf,
        MPI_Comm comm);

/**
 * Initialize an array given a starting seed.
 *
 * This function is used by both the client and the server to generate
 * initialize values in a buffer.
 */
extern void xfer_init_data_array(
        const unsigned int seed,
        data_array_t *array);

/**
 * Compare two arrays.  If they are equal, return 0.
 */
extern int xfer_compare_data_arrays(
        const data_array_t *arr1,
        const data_array_t *arr2);


extern int xfer_validate_array(
        const int seed,
        data_array_t *array);

extern int xfer_block_partition(
        const int num_bins,
        const int num_vals,
        const int val,
        int *bin,
        int *rank);


#ifdef __cplusplus
}
#endif

#endif /* XFER_UTIL_H_ */
