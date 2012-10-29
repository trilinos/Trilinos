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
 * aggregation.cpp
 *
 *  Created on: Mar 9, 2009
 *      Author: thkorde
 */

#include <assert.h>
#include <memory.h>

#include <Trios_nssi_server.h>  /* for nssi service functions */

#include "pnetcdf.h"
#include <mpi.h>
#include <algorithm>
#include <map>
#include <list>

using namespace std;

#include "netcdf_args.h"
#include "netcdf_debug.h"  /* netcdf_debug_level */

#include "aggregation.h"

typedef struct {
    int ahead_count;
    int behind_count;
    int same_count;
    int no_match_count;
} chunk_location_count_t;

struct participant_t {
    NNTI_peer_t p;
};
typedef struct participant_t participant_t;

typedef list<participant_t *> participants_t;
typedef list<participant_t *>::iterator participants_iterator_t;

typedef list<aggregation_chunk_details_t *> chunk_details_t;
typedef list<aggregation_chunk_details_t *>::iterator chunk_details_iterator_t;

typedef struct {
    aggregation_chunk_details_t *details;
    chunk_details_t              component_chunks;
} aggregation_chunk_t;

typedef list<aggregation_chunk_t *> chunks_t;
typedef list<aggregation_chunk_t *>::iterator chunks_iterator_t;


typedef struct {
    int varid;
    chunks_t *chunks;
    int       chunks_received;
} per_var_details_t;

typedef map<int, per_var_details_t *> nc_var_map_t;
typedef map<int, per_var_details_t *>::iterator nc_var_map_iterator_t;
typedef pair<int, per_var_details_t *> nc_var_map_pair_t;

typedef struct {
    int ncid;
    participants_t *participants;
    int             num_participants; /* will be init to 0.  if > 0, then someone forced a specific number of participants. */
    nc_var_map_t    nc_vars;
    write_type      type; /* direct, aggregate independent, aggregate collective */
} file_details_t;

static map<int, file_details_t *> open_file_map;
typedef map<int, file_details_t *>::iterator open_file_map_iterator_t;
typedef pair<int, file_details_t *> open_file_map_pair_t;

bool compare_chunks_for_aggregation(const aggregation_chunk_t* c1, const aggregation_chunk_t* c2)
{
    aggregation_chunk_details_t *details1=c1->details;
    aggregation_chunk_details_t *details2=c2->details;

    for (int i=0;i<details1->ndims;i++) {
        if (details1->count[i] < details2->count[i]) {
            return true;
        } else if (details1->count[i] > details2->count[i]) {
            return false;
        }
    }
    for (int i=0;i<details1->ndims;i++) {
        if (details1->start[i] < details2->start[i]) {
            return true;
        } else if (details1->start[i] > details2->start[i]) {
            return false;
        }
    }
    return false;
}

bool compare_chunks_for_caching(const aggregation_chunk_t* c1, const aggregation_chunk_t* c2)
{
    aggregation_chunk_details_t *details1=c1->details;
    aggregation_chunk_details_t *details2=c2->details;

    for (int i=0;i<details1->ndims;i++) {
        if (details1->start[i] < details2->start[i]) {
            return true;
        } else if (details1->start[i] > details2->start[i]) {
            return false;
        }
    }
    return false;
}

file_details_t *new_open_file(const int ncid)
{
    file_details_t *details=NULL;

    details = new file_details_t;
    details->ncid = ncid;
    details->participants = new participants_t;
    details->num_participants=0;

    return details;
}

int use_aggregation(const int ncid)
{
    file_details_t *details=NULL;

    details = open_file_map[ncid];
    if ((details->type == WRITE_AGGREGATE_INDEPENDENT) ||
        (details->type == WRITE_AGGREGATE_COLLECTIVE)) {
        return TRUE;
    }

    return FALSE;
}

int use_caching(const int ncid)
{
    file_details_t *details=NULL;

    details = open_file_map[ncid];
    if ((details->type == WRITE_CACHING_INDEPENDENT) ||
        (details->type == WRITE_CACHING_COLLECTIVE)) {
        return TRUE;
    }

    return FALSE;
}

int use_collective(const int ncid)
{
    file_details_t *details=NULL;

    details = open_file_map[ncid];
    if ((details->type == WRITE_AGGREGATE_COLLECTIVE) ||
        (details->type == WRITE_CACHING_COLLECTIVE)) {
        return TRUE;
    }

    return FALSE;
}

int use_independent(const int ncid)
{
    file_details_t *details=NULL;

    details = open_file_map[ncid];
    if ((details->type == WRITE_AGGREGATE_INDEPENDENT) ||
        (details->type == WRITE_CACHING_INDEPENDENT)) {
        return TRUE;
    }

    return FALSE;
}

int use_direct(const int ncid)
{
    file_details_t *details=NULL;

    details = open_file_map[ncid];
    if (details->type == WRITE_DIRECT) {
        return TRUE;
    }

    return FALSE;
}

void add_participant_for_file(const int ncid,
                              const NNTI_peer_t *caller,
                              const write_type write_type)
{
    file_details_t *details    =NULL;
    participant_t  *participant=NULL;

    log_debug(netcdf_debug_level, "adding participant for ncid=%d", ncid);
    details = open_file_map[ncid];
    if (details == NULL) {
        log_debug(netcdf_debug_level, "creating new file_details_t for ncid=%d", ncid);
        details=new_open_file(ncid);
        log_debug(netcdf_debug_level, "new file_details_t=%p", details);
        open_file_map[ncid]=details;
    }

    participant=new participant_t();
    participant->p=*caller;
    details->participants->push_back(participant);

    details->type = write_type;
}

void set_participant_count(const int ncid,
                           const NNTI_peer_t *caller,
                           const write_type write_type,
                           const int num_participants)
{
    file_details_t *details    =NULL;
    participant_t  *participant=NULL;

    log_debug(netcdf_debug_level, "setting participant count for ncid=%d", ncid);
    details = open_file_map[ncid];
    if (details == NULL) {
        log_debug(netcdf_debug_level, "creating new file_details_t for ncid=%d", ncid);
        details=new_open_file(ncid);
        log_debug(netcdf_debug_level, "new file_details_t=%p", details);
        open_file_map[ncid]=details;
    }

    participant=new participant_t();
    participant->p=*caller;
    details->participants->push_back(participant);
    details->type = write_type;
    details->num_participants=num_participants;
}

void remove_participant_for_file(const int ncid, const NNTI_peer_t *caller)
{
    file_details_t *details    =NULL;
    participant_t  *participant=NULL;
    participants_iterator_t iter;

    log_level debug_level = netcdf_debug_level;

    details = open_file_map[ncid];
    if (details == NULL) {
        log_error(debug_level, "remove failed - no participants");
        goto out;
    }

    iter = details->participants->begin();
    for (;iter != details->participants->end(); ++iter) {
        participant = *iter;
//        if (participant->p == *caller) {
//            details->participants->remove(participant);
//            delete(participant);
//            break;
//        }
    }

out:
    return;
}

void add_participant_chunk(const NNTI_peer_t *caller,
                           aggregation_chunk_details_t *chunk_details)
{
    int rc=NSSI_OK;

    file_details_t  *file_details=NULL;

    log_level debug_level = netcdf_debug_level;

    log_debug(netcdf_debug_level, "adding chunk: ncid(%d) varid(%d)", chunk_details->ncid, chunk_details->varid);

    chunk_details->datatype_size=0;
    MPI_Type_size(chunk_details->datatype, &chunk_details->datatype_size);
    if ((chunk_details->len/chunk_details->num_elements) != chunk_details->datatype_size) {
        log_warn(netcdf_debug_level, "datatype size conflict: (%d/%d)==%d is not equal to %d",
                 chunk_details->len, chunk_details->num_elements, chunk_details->len/chunk_details->num_elements, chunk_details->datatype_size);
    }

    file_details = open_file_map[chunk_details->ncid];
    if (file_details == NULL) {
        log_error(netcdf_debug_level, "failed to add chunk.  cannot aggregate.");
        return;
    }
    per_var_details_t *var_details = file_details->nc_vars[chunk_details->varid];
    if (var_details == NULL) {
        var_details=new per_var_details_t;
        var_details->varid = chunk_details->varid;
        var_details->chunks = new chunks_t;
        var_details->chunks_received=0;
        file_details->nc_vars[chunk_details->varid]=var_details;
    }
    aggregation_chunk_t *chunk=new aggregation_chunk_t;
    chunk->details = chunk_details;
    var_details->chunks->push_back(chunk);
    var_details->chunks_received++;

    return;
}

void destroy_chunk(aggregation_chunk_details_t *details)
{
    free(details->start);
    free(details->count);
    free(details->stride);
    log_debug(netcdf_debug_level, "freeing details->buf(%p)", details->buf);
    free(details->buf);
    delete details;
}

void cleanup_aggregation_chunks(const int ncid)
{
    file_details_t  *details=NULL;
    nc_var_map_iterator_t var_iter;
    per_var_details_t *var_details=NULL;

    log_debug(netcdf_debug_level, "entered");
    log_debug(netcdf_debug_level, "cleaning up - ncid(%d)", ncid);

    details = open_file_map[ncid];
    if (details == NULL) {
        return;
    }
    var_iter = details->nc_vars.begin();
    for (; var_iter != details->nc_vars.end(); ++var_iter) {
        var_details = var_iter->second;
        if (!var_details) {
            continue;
        }
        cleanup_aggregation_chunks(ncid, var_details->varid);
    }
}

void cleanup_aggregation_chunks(const int ncid, const int varid)
{
    file_details_t  *details=NULL;
    aggregation_chunk_t *chunk=NULL;
    chunks_iterator_t chunks_iter;
    chunk_details_iterator_t component_iter;
    nc_var_map_iterator_t vars_iter;

    log_debug(netcdf_debug_level, "cleaning up - ncid(%d) varid(%d)", ncid, varid);

    // for each variable, iterate over the chunks and destroy them

    details = open_file_map[ncid];

    per_var_details_t *var_details = details->nc_vars[varid];
    if (var_details != NULL) {
        chunks_iter = var_details->chunks->begin();
        for (;chunks_iter != var_details->chunks->end(); ++chunks_iter) {
            chunk = *chunks_iter;
            component_iter = chunk->component_chunks.begin();
            for (;component_iter != chunk->component_chunks.end(); ++component_iter) {
                log_debug(netcdf_debug_level, "cleanup - destroying component");
                destroy_chunk(*component_iter);
            }
            chunk->component_chunks.clear();
            log_debug(netcdf_debug_level, "cleanup - destroying details");
            destroy_chunk(chunk->details);
            delete chunk;
        }
        var_details->chunks->clear();
        var_details->chunks_received=0;
    }
}

static void recursive_print_chunk(aggregation_chunk_details_t *details, int offset, int *index, int current_dim)
{
    int my_offset=0;
    char tmp_str[20];
    char out_str[1024];
    int remaining=1023;

    if (current_dim < details->ndims-1) {
        for (int i=0;i<details->count[current_dim];i++) {
            my_offset = index[current_dim];
            for (int i=current_dim+1;i<details->ndims;i++) {
                my_offset *= details->count[i];
            }

            index[current_dim+1]=0;
            recursive_print_chunk(details, offset+my_offset, index, current_dim+1);
            index[current_dim] += details->datatype_size;
        }
        log_debug(netcdf_debug_level, "-----------------------------");
    } else {
        if (details->buf == NULL) {
            log_debug(netcdf_debug_level, "details->buf == NULL");
        } else {
            out_str[0]='\0';
            for (int i=0;i<details->count[current_dim];i++) {
                my_offset = offset+index[current_dim];

//                if (i==0) log_debug(netcdf_debug_level, "[%d][%d][%d] (my_offset==%d)", index[0], index[1], index[2], my_offset);
                if (details->datatype == MPI_BYTE || details->datatype == MPI_CHAR) {
                    sprintf(tmp_str, "%c, ", *(char *)(((char *)details->buf) + my_offset));
                }
                else if (details->datatype == MPI_SHORT) {
                    sprintf(tmp_str, "%hx, ", *(short *)(((char *)details->buf) + my_offset));
                }
                else if (details->datatype == MPI_INT) {
                    sprintf(tmp_str, "%x, ", *(int *)(((char *)details->buf) + my_offset));
                }
                else if (details->datatype == MPI_FLOAT) {
                    sprintf(tmp_str, "%f, ", *(float *)(((char *)details->buf) + my_offset));
                }
                else if (details->datatype == MPI_DOUBLE) {
                    sprintf(tmp_str, "%f, ", *(double *)(((char *)details->buf) + my_offset));
                }
                strncat(out_str, tmp_str, remaining);
                remaining -= strlen(out_str);

                index[current_dim] += details->datatype_size;
            }
//            log_debug(netcdf_debug_level, "[%d][%d][%d] (my_offset==%d)", index[0], index[1], index[2], my_offset);
            log_debug(netcdf_debug_level, "%s", out_str);
        }
    }
}

void print_chunk(aggregation_chunk_details_t *details)
{
    int *index=(int *)calloc(details->ndims, sizeof(int));
    int offset=0;
    int current_dim=0;
    char tmp_str[20];
    char out_str[1024];
    int remaining=1023;

    log_debug(netcdf_debug_level, "+++++++++++++++++++++++++++++");

    log_debug(netcdf_debug_level, "ncid==%d", details->ncid);
    log_debug(netcdf_debug_level, "varid==%d", details->varid);
    log_debug(netcdf_debug_level, "ndims==%d", details->ndims);
    log_debug(netcdf_debug_level, "len==%ld", details->len);
    log_debug(netcdf_debug_level, "num_elements==%d", details->num_elements);
    out_str[0]='\0';
    remaining=1023;
    for (int i=0;(i<details->ndims) && (remaining>0);i++) {
        sprintf(tmp_str, "%lld,", details->start[i]);
        strncat(out_str, tmp_str, remaining);
        remaining -= strlen(tmp_str);
    }
    log_debug(netcdf_debug_level, "start[]==%s", out_str);
    out_str[0]='\0';
    remaining=1023;
    for (int i=0;(i<details->ndims) && (remaining>0);i++) {
        sprintf(tmp_str, "%lld,", details->count[i]);
        strncat(out_str, tmp_str, remaining);
        remaining -= strlen(tmp_str);
    }
    log_debug(netcdf_debug_level, "count[]==%s", out_str);
    out_str[0]='\0';
    remaining=1023;
    for (int i=0;(i<details->ndims) && (remaining>0);i++) {
        sprintf(tmp_str, "%lld,", details->stride[i]);
        strncat(out_str, tmp_str, remaining);
        remaining -= strlen(tmp_str);
    }
    log_debug(netcdf_debug_level, "stride[]==%s", out_str);

//    recursive_print_chunk(details, offset, index, current_dim);
    log_debug(netcdf_debug_level, "+++++++++++++++++++++++++++++");

    free(index);
}

void print_chunk(aggregation_chunk_t *c)
{
    if (c->details == NULL) {
        log_debug(netcdf_debug_level, "chunk has no details.  perhaps it was aggregated into another chunk.");
        return;
    }
    print_chunk(c->details);
}

static void recursive_copy_chunk(aggregation_chunk_details_t *src,
                                 aggregation_chunk_details_t *dst,
                                 long src_offset,
                                 long dst_offset,
                                 int *src_index,
                                 int *dst_index,
                                 int current_dim)
{
    long my_src_offset=0;
    long my_dst_offset=0;

    if (current_dim < src->ndims-1) {
        for (int i=0;i<src->count[current_dim];i++) {
            my_src_offset = src_index[current_dim];
            my_dst_offset = dst_index[current_dim];
//            log_debug(netcdf_debug_level, "join_offset(%d) start_diff[%d](%d)",
//                    join_offset, current_dim, src->start[current_dim] - dst->start[current_dim]);
            my_dst_offset += ((src->start[current_dim] - dst->start[current_dim]) * src->datatype_size);
            for (int j=current_dim+1;j<src->ndims;j++) {
                my_src_offset *= src->count[j];
                my_dst_offset *= dst->count[j];
            }

            src_index[current_dim+1]=0;
            dst_index[current_dim+1]=0;
            recursive_copy_chunk(src, dst, src_offset+my_src_offset, dst_offset+my_dst_offset,
                                 src_index, dst_index, current_dim+1);
            src_index[current_dim] += src->datatype_size;
            dst_index[current_dim] += dst->datatype_size;
        }
    } else {
        dst_offset += ((src->start[current_dim] - dst->start[current_dim]) * src->datatype_size);
        memcpy(((char *)dst->buf) + dst_offset,
               ((char *)src->buf) + src_offset,
               src->count[current_dim]*src->datatype_size);
    }
}

static void recursive_aggregate_chunks(aggregation_chunk_t *src1,
                                       aggregation_chunk_t *src2,
                                       aggregation_chunk_t *dst)
{
    int *src_index=(int *)calloc(src1->details->ndims, sizeof(int));
    int *dst_index=(int *)calloc(dst->details->ndims, sizeof(int));
    int src_offset=0;
    long dst_offset=0;
    long current_dim=0;

    memset(src_index, 0, src1->details->ndims*sizeof(int));
    memset(dst_index, 0, dst->details->ndims*sizeof(int));
    recursive_copy_chunk(src1->details, dst->details, src_offset, dst_offset, src_index, dst_index, current_dim);
    memset(src_index, 0, src2->details->ndims*sizeof(int));
    memset(dst_index, 0, dst->details->ndims*sizeof(int));
    recursive_copy_chunk(src2->details, dst->details, src_offset, dst_offset, src_index, dst_index, current_dim);

    free(src_index);
    free(dst_index);
}

static void copy_chunk(aggregation_chunk_details_t *src,
                       aggregation_chunk_details_t *dst)
{
    int *src_index=(int *)calloc(src->ndims, sizeof(int));
    int *dst_index=(int *)calloc(dst->ndims, sizeof(int));
    int src_offset=0;
    long dst_offset=0;
    long current_dim=0;

    memset(src_index, 0, src->ndims*sizeof(int));
    memset(dst_index, 0, dst->ndims*sizeof(int));
    recursive_copy_chunk(src, dst, src_offset, dst_offset, src_index, dst_index, current_dim);

    free(src_index);
    free(dst_index);
}

aggregation_chunk_t *aggregate_chunks(aggregation_chunk_t *c1,
                                      aggregation_chunk_t *c2,
                                      int join_dim)
{
    aggregation_chunk_t *out=new aggregation_chunk_t;

    int src_buf_size=0;

    log_debug(netcdf_debug_level, "entered");

    assert(c1->details->ndims == c2->details->ndims);
    assert(out != NULL);

    out->details = new aggregation_chunk_details_t;

    out->details->ncid          = c1->details->ncid;
    out->details->varid         = c1->details->varid;
    out->details->ndims         = c1->details->ndims;
//    out->details->buf           = calloc(c1->details->len+c2->details->len, c1->details->datatype_size);
    out->details->buf           = NULL;
    out->details->atype         = c1->details->atype;
    out->details->len           = c1->details->len+c2->details->len;
    out->details->datatype      = c1->details->datatype;
    out->details->num_elements  = c1->details->num_elements+c2->details->num_elements;
    out->details->datatype_size = c1->details->datatype_size;
    out->details->start         = (MPI_Offset *)calloc(c1->details->ndims, sizeof(MPI_Offset));
    out->details->count         = (MPI_Offset *)calloc(c1->details->ndims, sizeof(MPI_Offset));
    out->details->stride        = (MPI_Offset *)calloc(c1->details->ndims, sizeof(MPI_Offset));

    memcpy(out->details->start, c1->details->start, c1->details->ndims*sizeof(MPI_Offset));
    memcpy(out->details->count, c1->details->count, c1->details->ndims*sizeof(MPI_Offset));
    out->details->count[join_dim] += c2->details->count[join_dim];
    memcpy(out->details->stride, c1->details->stride, c1->details->ndims*sizeof(MPI_Offset));

//    recursive_aggregate_chunks(c1, c2, out);

    if (c1->component_chunks.size() > 0) {
        out->component_chunks.merge(c1->component_chunks);
        c1->component_chunks.clear();
        destroy_chunk(c1->details);
    } else {
        out->component_chunks.push_back(c1->details);
    }
    c1->details = NULL;
    if (c2->component_chunks.size() > 0) {
        out->component_chunks.merge(c2->component_chunks);
        c2->component_chunks.clear();
        destroy_chunk(c2->details);
    } else {
        out->component_chunks.push_back(c2->details);
    }
    c2->details = NULL;

    assert(out != NULL);

    log_debug(netcdf_debug_level, "finished");

    return(out);
}

/*
 * Aggregate a particular variable in the file.
 *
 * Aggregation rules:
 *  - dimension count must be equal
 *  - strides must be equal
 *  - counts on matching faces must be equal
 *  -
 *
 */
int try_aggregation(const int ncid, const int varid)
{
    int aggregation_success=FALSE;

    file_details_t  *file_details=NULL;
    per_var_details_t *var_details=NULL;
    aggregation_chunk_t *base_chunk=NULL;
    aggregation_chunk_t *candidate_chunk=NULL;
    aggregation_chunk_t *new_chunk=NULL;
    chunks_iterator_t base_iter, candidate_iter;
    int *start_diff;
    chunk_location_count_t chunk_location_count;
    int dim_with_movement=-1;

    chunks_t agg_chunks;

    int failed=FALSE;

    file_details = open_file_map[ncid];
    if (file_details == NULL) {
        return(aggregation_success);
    }
    var_details = file_details->nc_vars[varid];
    if (var_details == NULL) {
        return(aggregation_success);
    }
    if (var_details->chunks->size() < 2) {
        log_debug(netcdf_debug_level, "returning with chunk count(%d)", var_details->chunks->size());
        return(aggregation_success);
    }
    log_debug(netcdf_debug_level, "chunk count(%d)", var_details->chunks->size());


    log_debug(netcdf_debug_level, "trying aggregation - ncid(%d) varid(%d)", ncid, varid);

    var_details->chunks->sort(compare_chunks_for_aggregation);

    log_level old=netcdf_debug_level;
//    netcdf_debug_level=LOG_ALL;
//    log_debug(netcdf_debug_level, "*****************");
//    log_debug(netcdf_debug_level, "start aggregation (begin list)");
//    log_debug(netcdf_debug_level, "*****************");
//    int chunk_count;
//    aggregation_chunk_t **chunks = get_chunks(ncid, varid, &chunk_count);
//    for (int i=0;i<chunk_count;i++) {
//        print_chunk(chunks[i]);
//    }
//    free(chunks);
//    log_debug(netcdf_debug_level, "*****************");
//    log_debug(netcdf_debug_level, "start aggregation (end list)");
//    log_debug(netcdf_debug_level, "*****************");
//    netcdf_debug_level=old;

    int success_this_pass=TRUE;
    while (success_this_pass==TRUE) {
        success_this_pass=FALSE;

//        log_debug(LOG_ALL, "top: while loop");


        base_iter = var_details->chunks->begin();
        base_chunk = *base_iter;
        start_diff=new int[base_chunk->details->ndims];
        for (;base_iter != var_details->chunks->end(); ++base_iter) {
//            log_debug(LOG_ALL, "top: base_iter loop");

            base_chunk = *base_iter;

            //if (base_chunk != NULL)      print_chunk(base_chunk);

            // look for a chunk that can be aggregated to the base chunk
            candidate_iter = base_iter;
            candidate_iter++;
            for (;candidate_iter != var_details->chunks->end(); ++candidate_iter) {
//                log_debug(LOG_ALL, "top: candidate_iter loop");

                candidate_chunk = *candidate_iter;

                //if (candidate_chunk != NULL) print_chunk(candidate_chunk);

                failed=FALSE;

                if (base_chunk->details->ndims != candidate_chunk->details->ndims) {
                    continue;
                }
//                if (candidate_chunk->details->start[0] != base_chunk->details->start[0]) {
//                    continue;
//                }
                for (int i=0; i<base_chunk->details->ndims; i++) {
                    if (base_chunk->details->stride[i] != candidate_chunk->details->stride[i]) {
                        failed=TRUE;
                        break;
                    }

                    start_diff[i] = candidate_chunk->details->start[i] - base_chunk->details->start[i];
                }
                if (failed) continue;

                chunk_location_count.ahead_count=0;
                chunk_location_count.behind_count=0;
                chunk_location_count.same_count=0;
                chunk_location_count.no_match_count=0;
                int agg_dims=base_chunk->details->ndims; /* the number of dimensions to aggregate */
                int first_agg_dim=0; /* first dimensions to aggregate */
                for (int i=first_agg_dim; i<agg_dims; i++) {
                    if ((start_diff[i] < 0) && (-start_diff[i] == candidate_chunk->details->count[i])) {
                        // the candidate is "behind/below" and touching the base chunk in this dimension
                        chunk_location_count.behind_count++;
                        dim_with_movement=i;
                    } else if ((start_diff[i] > 0) && (start_diff[i] == base_chunk->details->count[i])) {
                        // the candidate is "ahead of/above" and touching the base chunk in this dimension
                        chunk_location_count.ahead_count++;
                        dim_with_movement=i;
                    } else if (start_diff[i] == 0) {
                        // the candidate is "equal to" the base chunk in this dimension
                        chunk_location_count.same_count++;
                    } else {
                        // the candidate and the base chunk don't match in this dimension
                        chunk_location_count.no_match_count++;
                    }
                }

#ifdef DEBUG
                /*
                 * These tests can be interesting, but are not required to get the job done.
                 */
                if (chunk_location_count.no_match_count > 0) {
                    // no matching face found.  can't aggregate.
                    continue;
                }

                if (chunk_location_count.same_count == base_chunk->ndims) {
                    // base and candidate have same start.  bad?  can't aggregate.
                    continue;
                }

                if (chunk_location_count.ahead_count > 1) {
                    // movement in more than one direction
                    continue;
                }
                if (chunk_location_count.behind_count > 1) {
                    // movement in more than one direction
                    continue;
                }
                if ((chunk_location_count.ahead_count > 0)  &&
                        (chunk_location_count.behind_count > 0)) {
                    // movement in more than one direction
                    continue;
                }

                if ((chunk_location_count.ahead_count == 0)  &&
                        (chunk_location_count.behind_count == 0)) {
                    // possible movement, but the chunks don't touch
                    continue;
                }
#endif

                // check that the matching faces have the same dimensions
                for (int i=0; i<base_chunk->details->ndims; i++) {
                    if ((i != dim_with_movement) &&
                        (base_chunk->details->count[i] != candidate_chunk->details->count[i])) {
                        failed=TRUE;
                        break;
                    }
                }
                if (failed) continue;

                /*
                 * Do NOT uncomment these print_chunk() lines in production code.
                 * They are *very* slow even if the debug level is set low and
                 * nothing is being logged.
                 */
//                netcdf_debug_level=LOG_ALL;
//                log_debug(netcdf_debug_level, "*****************");
//                log_debug(netcdf_debug_level, "base chunk");
//                log_debug(netcdf_debug_level, "*****************");
//                if (base_chunk != NULL)      print_chunk(base_chunk);
//                log_debug(netcdf_debug_level, "*****************");
//                log_debug(netcdf_debug_level, "candidate chunk");
//                log_debug(netcdf_debug_level, "*****************");
//                if (candidate_chunk != NULL) print_chunk(candidate_chunk);
//                netcdf_debug_level=old;

                if ((chunk_location_count.ahead_count == 1)  &&
                        (chunk_location_count.behind_count == 0) &&
                        (chunk_location_count.same_count == agg_dims-1)) {
                    // aggregation is base + candidate
                    new_chunk = aggregate_chunks(base_chunk, candidate_chunk, dim_with_movement);
                } else if ((chunk_location_count.ahead_count == 0)  &&
                        (chunk_location_count.behind_count == 1) &&
                        (chunk_location_count.same_count == agg_dims-1)) {
                    // aggregation is candidate + base
                    new_chunk = aggregate_chunks(candidate_chunk, base_chunk, dim_with_movement);
                } else {
                    // chunks aren't aligned
                    //printf("**********\nchunks are not aligned\n**********\n");
                    continue;
                }

                assert(new_chunk != NULL);

                /*
                 * Do NOT uncomment these print_chunk() lines in production code.
                 * They are *very* slow even if the debug level is set low and
                 * nothing is being logged.
                 */
//                netcdf_debug_level=LOG_ALL;
//                log_debug(netcdf_debug_level, "*****************");
//                log_debug(netcdf_debug_level, "new chunk");
//                log_debug(netcdf_debug_level, "*****************");
//                if (new_chunk != NULL)       print_chunk(new_chunk);
//                netcdf_debug_level=old;

                var_details->chunks->remove(base_chunk);
                var_details->chunks->remove(candidate_chunk);
                delete base_chunk;
                delete candidate_chunk;

                agg_chunks.push_back(new_chunk);

                aggregation_success = TRUE;
                success_this_pass = TRUE;

                break;
            }
            if (success_this_pass == TRUE) break;
        }
        chunks_iterator_t agg_iter = agg_chunks.begin();
        for (;agg_iter != agg_chunks.end();agg_iter++) {
            var_details->chunks->push_back(*agg_iter);
        }
        agg_chunks.clear();

        delete[] start_diff;
    }

//    netcdf_debug_level=LOG_ALL;
//    log_debug(netcdf_debug_level, "*****************");
//    log_debug(netcdf_debug_level, "end aggregation (begin list)");
//    log_debug(netcdf_debug_level, "*****************");
//    chunks = get_chunks(ncid, varid, &chunk_count);
//    for (int i=0;i<chunk_count;i++) {
//        print_chunk(chunks[i]);
//    }
//    free(chunks);
//    log_debug(netcdf_debug_level, "*****************");
//    log_debug(netcdf_debug_level, "end aggregation (end list)");
//    log_debug(netcdf_debug_level, "*****************");
//    netcdf_debug_level=old;

//    netcdf_debug_level=LOG_ALL;
    chunks_iterator_t dst_iter=var_details->chunks->begin();
    for(;dst_iter != var_details->chunks->end();dst_iter++) {
        chunk_details_iterator_t component_iter=(*dst_iter)->component_chunks.begin();
        if (((*dst_iter)->details->buf == NULL) && ((*dst_iter)->details->len > 0)) {
            (*dst_iter)->details->buf = (char *)malloc((*dst_iter)->details->len);
            log_debug(netcdf_debug_level, "allocated dst_iter->details->buf(%p), len(%ld)",
                    (*dst_iter)->details->buf,
                    (*dst_iter)->details->len);
        } else {
            log_debug(netcdf_debug_level, "did not allocate dst_iter->details->buf(%p)", (*dst_iter)->details->buf);
        }
        for(;component_iter != (*dst_iter)->component_chunks.end();component_iter++) {
            log_debug(netcdf_debug_level, "copying component");
            copy_chunk(*component_iter, (*dst_iter)->details);
            log_debug(netcdf_debug_level, "destroying component");
            destroy_chunk(*component_iter);
        }
        (*dst_iter)->component_chunks.clear();
    }
//    netcdf_debug_level=old;

//    netcdf_debug_level=LOG_ALL;
    log_debug(netcdf_debug_level, "*****************");
    log_debug(netcdf_debug_level, "chunks after aggregation");
    log_debug(netcdf_debug_level, "*****************");
    base_iter = var_details->chunks->begin();
    for (;base_iter != var_details->chunks->end(); ++base_iter) {
        base_chunk = *base_iter;
        if (base_chunk != NULL)
            print_chunk(base_chunk);
    }
//    netcdf_debug_level=old;

    return(aggregation_success);
}

/*
 * Aggregate all variables in the file.
 *
 */
int try_aggregation(const int ncid)
{
    int aggregation_success=FALSE;

    file_details_t  *file_details=NULL;
    nc_var_map_iterator_t var_iter;
    per_var_details_t *var_details=NULL;

    file_details = open_file_map[ncid];
    if (file_details == NULL) {
        return(aggregation_success);
    }
    var_iter = file_details->nc_vars.begin();
    for (; var_iter != file_details->nc_vars.end(); var_iter++) {
        var_details = var_iter->second;
        if (!var_details) {
            continue;
        }
        while(try_aggregation(ncid, var_details->varid) == TRUE);
    }

    aggregation_success = TRUE;

    return(aggregation_success);
}

int aggregate_data_ready_to_write(const int ncid, const int varid)
{
//    file_details_t *details = open_file_map[ncid];
//    int chunks_needed=0;

//    if (details->num_participants > 0) {
//        chunks_needed = details->num_participants;
//    } else {
//        chunks_needed = details->participants->size();
//    }
//
//    if (details->nc_vars[varid]->chunks_received == chunks_needed) {
//        return TRUE;
//    }

    return FALSE;
}

int cache_data_ready_to_write(const int ncid, const int varid)
{
//    file_details_t *details = open_file_map[ncid];
//    int chunks_needed=0;

//    if (details->num_participants > 0) {
//        chunks_needed = details->num_participants;
//    } else {
//        chunks_needed = details->participants->size();
//    }
//
//    if (details->nc_vars[varid]->chunks_received == chunks_needed) {
//        return TRUE;
//    }

    return FALSE;
}

aggregation_chunk_details_t **get_chunks(const int ncid, const int varid, int *chunk_count)
{
    file_details_t *details=NULL;
    per_var_details_t *var_details=NULL;
    aggregation_chunk_details_t **chunks=NULL;
    chunks_iterator_t iter;

    log_debug(netcdf_debug_level, "entered");

    *chunk_count=0;

    details = open_file_map[ncid];
    if (details == NULL) {
        return(NULL);
    }
    var_details = details->nc_vars[varid];
    if (var_details == NULL) {
        return(NULL);
    }

    *chunk_count = details->nc_vars[varid]->chunks->size();

    log_debug(netcdf_debug_level, "found %d chunks to return", *chunk_count);

    if (*chunk_count == 0) {
        return(NULL);
    }
    chunks = (aggregation_chunk_details_t **)malloc(*chunk_count*sizeof(aggregation_chunk_details_t *));

    var_details->chunks->sort(compare_chunks_for_caching);

    iter = var_details->chunks->begin();
    for (int i=0;iter != var_details->chunks->end(); ++iter,i++) {
        chunks[i] = (*iter)->details;
//        print_chunk(chunks[i]);
    }

    log_debug(netcdf_debug_level, "finished");

    return(chunks);
}

aggregation_chunk_details_t **get_chunks(const int ncid, int *chunk_count)
{
    file_details_t  *details=NULL;
    nc_var_map_iterator_t var_iter;
    per_var_details_t *var_details=NULL;
    aggregation_chunk_details_t **chunks=NULL;
    chunks_iterator_t chunks_iter;

    log_debug(netcdf_debug_level, "entered");

    *chunk_count=0;

    details = open_file_map[ncid];
    if (details == NULL) {
        return(NULL);
    }
    var_iter = details->nc_vars.begin();
    for (; var_iter != details->nc_vars.end(); ++var_iter) {
        var_details = var_iter->second;
        if (!var_details) {
            continue;
        }
        *chunk_count += var_details->chunks->size();
    }

    log_debug(netcdf_debug_level, "found %d chunks to return", *chunk_count);

    if (*chunk_count == 0) {
        return(NULL);
    }
    chunks = (aggregation_chunk_details_t **)malloc(*chunk_count*sizeof(aggregation_chunk_details_t *));

    int i=0;
    var_iter = details->nc_vars.begin();
    for (; var_iter != details->nc_vars.end(); var_iter++) {
        var_details = var_iter->second;
        if (!var_details) {
            continue;
        }
        var_details->chunks->sort(compare_chunks_for_caching);
        chunks_iter = var_details->chunks->begin();
        for (;chunks_iter != var_details->chunks->end(); ++chunks_iter,i++) {
            chunks[i] = (*chunks_iter)->details;
//            print_chunk(chunks[i]);
        }
    }

    log_debug(netcdf_debug_level, "finished");

    return(chunks);
}
