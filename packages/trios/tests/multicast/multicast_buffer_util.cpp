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
 * multicast_init_buffer.cpp
 *
 *  Created on: Nov 4, 2011
 *      Author: raoldfi
 */


#include "multicast_service_args.h"
#include "multicast_debug.h"
#include <stdlib.h>


/**
 * Initialize an array given a starting seed.
 *
 * This function is used by both the client and the server to generate
 * initialize values in a buffer.
 */
void multicast_init_data_array(const unsigned int seed, data_array_t *array)
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
int multicast_compare_data_arrays(const data_array_t *arr1, const data_array_t *arr2)
{
    log_level debug_level = multicast_debug_level;

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


int multicast_validate_array(const int seed,  data_array_t *array)
{
    int rc = 0;
    log_level debug_level = multicast_debug_level;

    /* Validate the array that was sent through the args */
    data_array_t tmp_array;

    tmp_array.data_array_t_len = array->data_array_t_len;
    tmp_array.data_array_t_val = new data_t[array->data_array_t_len];

    multicast_init_data_array(seed, &tmp_array);
    rc = multicast_compare_data_arrays(&tmp_array, array);

    if (rc != 0) {
        log_warn(debug_level, "Unable to validate array");
    }

    delete [] tmp_array.data_array_t_val;
    return rc;
}
