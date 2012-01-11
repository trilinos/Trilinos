/*
 * xfer_init_buffer.cpp
 *
 *  Created on: Nov 4, 2011
 *      Author: raoldfi
 */


#include "xfer_service_args.h"
#include <stdlib.h>


/**
 * Initialize an array given a starting seed.
 *
 * This function is used by both the client and the server to generate
 * initialize values in a buffer.
 */
void xfer_init_data_array(const unsigned int seed, data_array_t *array)
{
    char state[8];

    int len = array->data_array_t_len;
    data_t *buf = array->data_array_t_val;

    /* initialize the random seed */
    initstate(seed, state, 8);

    for (int i=0; i<len; i++) {

        long rand_val = random();

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
    if (arr1->data_array_t_len != arr2->data_array_t_len)
        return -1;

    for (int i=0; i<arr1->data_array_t_len; i++) {
        if (arr1->data_array_t_val[i].int_val != arr2->data_array_t_val[i].int_val)
            return -1;
        if (arr1->data_array_t_val[i].float_val != arr2->data_array_t_val[i].float_val)
            return -1;
        if (arr1->data_array_t_val[i].double_val != arr2->data_array_t_val[i].double_val)
            return -1;
    }

    return 0;

}
