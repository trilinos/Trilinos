/*
 * xfer_util.cpp
 *
 *  Created on: Nov 4, 2011
 *      Author: raoldfi
 */


#include "Trios_nssi_types_xdr.h"
#include "xfer_service_args.h"
#include "xfer_debug.h"
#include "xfer_util.h"
#include "xfer_client.h"

#include <stdlib.h>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <cstdlib>



/**
 * Output the configuration file with this format
 * ------
 * line 1: ns  # number of servers
 * line 2: server_url_0
 * line 3: server_url_1
 * ...
 * line n: server_url_n-1
 */
int xfer_write_server_url_file(
        std::string url_fname,
        std::string my_url,
        MPI_Comm comm)
{
    int rc = 0;
    int rank, np;
    log_level debug_level = LOG_UNDEFINED;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);

    // gather the string lengths for the input strings
    int total_length=0;
    int max_length=0;
    int len = my_url.length();
    MPI_Reduce(&len, &total_length, 1, MPI_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&len, &max_length, 1, MPI_INT, MPI_MAX, 0, comm);

    // gather the lengths of the individual strings as offsets
    std::vector<int> length_vector(np, 0);
    MPI_Gather(&len, 1, MPI_INT, &length_vector[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

    // output the config file
    if (!url_fname.empty()) {

        // gather the urls from all the servers
        std::string urls(max_length*np,0);

        MPI_Gather(&my_url[0], len, MPI_CHAR, &urls[0], max_length, MPI_CHAR, 0, comm);

        // Write the urls to the config file
        if (rank == 0) {
            int i;

            std::ofstream urlfile (url_fname.c_str());
            if (urlfile.is_open()) {
                // write the number of servers as the first line
                urlfile << np << std::endl;
                size_t pos = 0;

                for (i=0; i<np; i++) {
                    std::string url = urls.substr(pos, length_vector[i]);
                    pos += length_vector[i];
                    log_debug(debug_level, "-- server %d : %s", i, url.c_str());
                    urlfile << url.c_str() << std::endl;
                }
            }
            else {
                log_error(debug_level, "Could not open server_url_file %s", url_fname.c_str());
                MPI_Abort(comm, 1);
            }

            // close the file
            urlfile.close();
        }


    }

    return rc;
}

/**
 * Get the list of server URLs from the file.
 */
void xfer_read_server_url_file(const char *path, std::vector<std::string> &urlbuf, MPI_Comm comm)
{
    log_level debug_level = xfer_debug_level;
    int rank, np;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);

    // variables for the serverURLs
    int numServers = 0;

    // Open the file and get all the server urls
    if (rank == 0) {

        // try to open the file
        std::ifstream urlFile(path);
        if (urlFile.is_open()) {
            if (urlFile.good()) {
                std::string line;
                std::getline(urlFile, line);
                numServers = std::atoi(line.c_str());
            }

            urlbuf.resize(numServers, std::string(NSSI_URL_LEN, '\0'));

            for (int i=0; i<numServers; i++) {
                if (urlFile.good()) {
                    std::getline(urlFile, urlbuf[i]);
                    log_debug(debug_level, "URL %d = %s", i, urlbuf[i].c_str());
                }
                else {
                    log_error(debug_level, "Unexpected EOF in %s", path);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
        }

        else {
            log_error(debug_level, "Could not open serverURL = %s", path);
            MPI_Abort(comm, 1);
        }
        urlFile.close();
    }

    // broadcast the number of servers
    MPI_Bcast(&numServers, 1, MPI_INT, 0, comm);

    if (rank != 0)
        urlbuf.resize(numServers, std::string(NSSI_URL_LEN, '\0'));  // allocate space for URLs


    // broadcast the list of URls to the other clients
    for (int i=0; i<numServers; i++) {
        MPI_Bcast(&urlbuf[i][0], NSSI_URL_LEN, MPI_CHAR, 0, comm);
    }
}

/**
 * Initialize an array given a starting seed.
 *
 * This function is used by both the client and the server to generate
 * initialize values in a buffer.
 */
void xfer_init_data_array(const unsigned int seed, data_array_t *array)
{
    //char state[8];

    int len = array->data_array_t_len;
    data_t *buf = array->data_array_t_val;

    /* initialize the random seed */
    //initstate(seed, state, 8);

    for (int i=0; i<len; i++) {

        //long rand_val = random();
        long rand_val = seed+1024+i;

        buf[i].int_val = (int)rand_val;
        buf[i].float_val = (float)rand_val;
        buf[i].double_val = (double)rand_val;
    }
}


/**
 * Compare two arrays.  If they are equal, return 0.
 */
int xfer_compare_data_arrays(const data_array_t *arr1, const data_array_t *arr2)
{
    log_level debug_level = xfer_debug_level;

    if (arr1->data_array_t_len != arr2->data_array_t_len) {
        log_error(debug_level, "arr1->len=%d, arr2->len=%d",
                arr1->data_array_t_len, arr2->data_array_t_len);
        return -1;
    }

    for (int i=0; i<(int)arr1->data_array_t_len; i++) {
        if (arr1->data_array_t_val[i].int_val != arr2->data_array_t_val[i].int_val) {
            log_error(debug_level, "val[%d].int_val=%d, val[%d].int_val=%d",
                    i,arr1->data_array_t_val[i].int_val,
                    i,arr2->data_array_t_val[i].int_val);
            return -1;
        }
        if (arr1->data_array_t_val[i].float_val != arr2->data_array_t_val[i].float_val) {
            log_error(debug_level, "val[%d].float_val=%f, val[%d].float_val=%f",
                                i,arr1->data_array_t_val[i].float_val,
                                i,arr2->data_array_t_val[i].float_val);
            return -1;
        }
        if (arr1->data_array_t_val[i].double_val != arr2->data_array_t_val[i].double_val) {
            log_error(debug_level, "val[%d].double_val=%g, val[%d].double_val=%g",
                                i,arr1->data_array_t_val[i].double_val,
                                i,arr2->data_array_t_val[i].double_val);
            return -1;
        }
    }

    return 0;

}


int xfer_validate_array(const int seed,  data_array_t *array)
{
    int rc = 0;
    log_level debug_level = xfer_debug_level;

    /* Validate the array that was sent through the args */
    data_array_t tmp_array;

    tmp_array.data_array_t_len = array->data_array_t_len;
    tmp_array.data_array_t_val = new data_t[array->data_array_t_len];

    xfer_init_data_array(seed, &tmp_array);
    rc = xfer_compare_data_arrays(&tmp_array, array);

    if (rc != 0) {
        log_warn(debug_level, "Unable to validate array");
    }

    delete [] tmp_array.data_array_t_val;
    return rc;
}


/**
 * Partition this val into a bin using a simple linear block-partitioning
 * scheme.
 *
 * Returns the bin assignment and the rank within that bin.
 *
 */
int xfer_block_partition(
        const int num_bins,
        const int num_vals,
        const int val,
        int *bin,
        int *rank)
{

    int per_bin = num_vals / num_bins;
    int num_large = num_vals % num_bins;  // number of bins with one extra

    int barrier = num_large * (per_bin + 1);  // total number of vals in large bins

    if (val < barrier) {
        *bin = val / (per_bin + 1);
        *rank = val % (per_bin + 1);
    }

    else {
        *bin = num_large + (val - barrier)/per_bin;
        *rank = (val - barrier) % per_bin;
    }

    //std::cout << "num_bins=" << num_bins << ", num_vals=" << num_vals;
    //std::cout << ": val=" << val << ", bin=" << bin << ", rank=" << rank << std::endl;

    return 0;
}
