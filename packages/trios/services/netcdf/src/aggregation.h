/*
 * aggregation.h
 *
 *  Created on: Mar 9, 2009
 *      Author: thkorde
 */

#ifndef AGGREGATION_H_
#define AGGREGATION_H_

#include "netcdf_args.h"


struct participant_t;
struct aggregation_chunk_details_t {
    participant_t *p;

    int ncid;
    int varid;
    int ndims;

    void *buf; /* the data */

    int    atype; /* netcdf type of data in buf*/
    size_t len;   /* length of buf in bytes */

    MPI_Datatype datatype; /* MPI type of data in buf */
    int num_elements;      /* number of datatype elements in buf (len/sizeof(datatype)) */
    int datatype_size;

    MPI_Offset *start;  /* starting corner (eg. 0,0,0 is the origin of a cube) */
    MPI_Offset *count;  /* num elements in each dimension (eg. 3,3,3 is a cube of size 3) */
    MPI_Offset *stride; /* sampling interval (eg. 2,2,2 is every other element in each dimension) */
};
typedef struct aggregation_chunk_details_t aggregation_chunk_details_t;

int use_aggregation(const int ncid);
int use_caching(const int ncid);
int use_collective(const int ncid);
int use_independent(const int ncid);
int use_direct(const int ncid);
void add_participant_for_file(const int ncid,
                              const NNTI_peer_t *caller,
                              const write_type write_type);
void set_participant_count(const int ncid,
                           const NNTI_peer_t *caller,
                           const write_type write_type,
                           const int num_participants);
void remove_participant_for_file(const int ncid, const NNTI_peer_t *caller);
void add_participant_chunk(const NNTI_peer_t *caller,
                           aggregation_chunk_details_t *chunk);
void cleanup_aggregation_chunks(const int ncid);
void cleanup_aggregation_chunks(const int ncid, const int varid);
int try_aggregation(const int ncid);
int try_aggregation(const int ncid, const int varid);
int aggregate_data_ready_to_write(const int ncid, const int varid);
int cache_data_ready_to_write(const int ncid, const int varid);
aggregation_chunk_details_t **get_chunks(const int ncid, int *chunk_count);
aggregation_chunk_details_t **get_chunks(const int ncid, const int varid, int *chunk_count);
void print_chunk(aggregation_chunk_details_t *c);



#endif /* AGGREGATION_H_ */
