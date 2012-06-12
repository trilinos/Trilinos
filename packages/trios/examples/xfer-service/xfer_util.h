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
