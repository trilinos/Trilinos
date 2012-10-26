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
#include <stdlib.h>
#include <string.h>

#include <netcdf.h>

#include <algorithm>
#include <map>
#include <list>
#include <vector>

using namespace std;

#include "create_subchunks.h"

#include "netcdf_debug.h"
#include "netcdf_config_parser.h"

#include <mpi.h>
#include "io_timer.h"


extern struct netcdf_config nc_cfg;


/*
 * if you want to be completely generic, MAX_SUBCHUNKING_DIMS should be
 * equal to NC_MAX_DIMS.  However, this causes a boat-load of memory to
 * allocated.  Then as you traverse the array, the memory must be mapped
 * into RAM which is SLOW.  Most codes we see use at most 4 dims.  Use 6
 * here just in case.
 */
//#define MAX_SUBCHUNKING_DIMS NC_MAX_DIMS
#define MAX_SUBCHUNKING_DIMS 6
struct start_offset_length_service {
    nc_size_t start[MAX_SUBCHUNKING_DIMS];  /* coordinates of a contiguous subchunk */
    nc_size_t count[MAX_SUBCHUNKING_DIMS];  /* dimensions of a contiguous subchunk */
    int       ndims;
    nc_size_t offset;
    nc_size_t length;
    int       datatype_size;
    int       service;
};
typedef struct start_offset_length_service start_offset_length_service_t;


/*
 * This is a special version of calc_offset_length() for the case where
 * there are 4 dimensions and highest_contig_dim==3.
 */
static void calc_offset_length_4dims(const size_t    *dimlens,
                                     const int        highest_contig_dim,
                                     start_offset_length_service_t *sols,
                                     int             *ol_index,
                                     const nc_size_t *start,
                                     const nc_size_t *count,
                                     const size_t     datatype_size)
{
    nc_size_t file_offset=0;
    nc_size_t lower_dim_products[4];


    lower_dim_products[3]=1;
    lower_dim_products[2]=dimlens[3];
    lower_dim_products[1]=dimlens[2]*dimlens[3];
    lower_dim_products[0]=dimlens[1]*dimlens[2]*dimlens[3];


    for (int i=0;i<count[0];i++) {
        for (int j=0;j<count[1];j++) {
            for (int k=0;k<count[2];k++) {
                // populate the sols for this subchunk
                sols[*ol_index].start[0]=start[0]+i;
                sols[*ol_index].start[1]=start[1]+j;
                sols[*ol_index].start[2]=start[2]+k;
                sols[*ol_index].start[3]=start[3];

                sols[*ol_index].count[0]=1;
                sols[*ol_index].count[1]=1;
                sols[*ol_index].count[2]=1;
                sols[*ol_index].count[3]=count[3];

                sols[*ol_index].ndims=4;
                sols[*ol_index].offset=(i*lower_dim_products[0]) +
                                       (j*lower_dim_products[1]) +
                                       (k*lower_dim_products[2]);
                sols[*ol_index].length=count[3]*datatype_size;
                sols[*ol_index].datatype_size=datatype_size;
                (*ol_index)++;
            }
        }
    }
//    log_debug(LOG_ALL, "ndims=%d; start[0,1,2,3]=%d,%d,%d,%d; count[0,1,2,3]=%d,%d,%d,%d; highest_contig_dim=%d, "
//            "current_dim=%d, file_offset=%d, current_start[0,1,2,3]=%d,%d,%d,%d; ol_index=%d",
//            ndims, start[0], start[1], start[2], start[3], count[0], count[1], count[2], count[3], highest_contig_dim,
//            current_dim, file_offset, current_start[0], current_start[1], current_start[2], current_start[3], *ol_index);
}

static void calc_offset_length(const int        ndims,
                               const size_t    *dimlens,
                               const int        highest_contig_dim,
                               const int        current_dim,
                               const nc_size_t  file_offset,
                               nc_size_t       *current_start,
                               start_offset_length_service_t *sols,
                               int             *ol_index,
                               const nc_size_t *start,
                               const nc_size_t *count,
                               const size_t     datatype_size)
{
//    log_debug(LOG_ALL, "ndims=%d; start[0,1,2,3]=%d,%d,%d,%d; count[0,1,2,3]=%d,%d,%d,%d; highest_contig_dim=%d, "
//            "current_dim=%d, file_offset=%d, current_start[0,1,2,3]=%d,%d,%d,%d; ol_index=%d",
//            ndims, start[0], start[1], start[2], start[3], count[0], count[1], count[2], count[3], highest_contig_dim,
//            current_dim, file_offset, current_start[0], current_start[1], current_start[2], current_start[3], *ol_index);
    if (current_dim < highest_contig_dim) {
        current_start[current_dim]=start[current_dim];
        for (int i=0;i<count[current_dim];i++) {
            nc_size_t my_offset_adder = current_start[current_dim];
            for (int j=current_dim+1;j<ndims;j++) {
                my_offset_adder *= dimlens[j];
            }
            my_offset_adder *= datatype_size;

            current_start[current_dim+1]=0;
            calc_offset_length(ndims, dimlens, highest_contig_dim, current_dim+1, file_offset+my_offset_adder,
                               current_start, sols, ol_index, start, count, datatype_size);
            current_start[current_dim]++;
        }
    } else {
        current_start[current_dim]=start[current_dim];
        nc_size_t my_offset_adder = start[current_dim];
        for (int j=current_dim+1;j<ndims;j++) {
            my_offset_adder *= dimlens[j];
        }
        my_offset_adder *= datatype_size;
        nc_size_t length=1;
        for (int j=current_dim;j<ndims;j++) {
            length *= count[j];
        }
        length *= datatype_size;

        // populate the sols for this subchunk
        memcpy(sols[*ol_index].start, current_start, ndims*sizeof(nc_size_t));
        for (int j=0;j<current_dim;j++) {
            sols[*ol_index].count[j]=1;
        }
        for (int j=current_dim;j<ndims;j++) {
            sols[*ol_index].count[j]=count[j];
        }
        sols[*ol_index].ndims=ndims;
        sols[*ol_index].offset=file_offset+my_offset_adder;
        sols[*ol_index].length=length;
        sols[*ol_index].datatype_size=datatype_size;
        (*ol_index)++;


//        log_debug(netcdf_debug_level,"start[");
//        for (int k=0;k<ndims;k++) {
//            log_debug(netcdf_debug_level,"%03lu,", current_start[k]);
//        }
//        log_debug(netcdf_debug_level,"] file_offset(%lu) length(%lu)", file_offset+my_offset_adder, length);
    }
}

static void assign_service(const nc_size_t dim_product,
                           const nc_size_t bytes_per_service,
                           start_offset_length_service_t *sols,
                           const int num_sols)
{
    for(int i=0;i<num_sols;i++) {
        sols[i].service = (sols[i].offset % dim_product)/bytes_per_service;
        log_debug(netcdf_debug_level, "sols[%d].offset(%lu) bytes_per_server(%lu) sols[%d].service(%d)",
                i, sols[i].offset, bytes_per_service, i, sols[i].service);
    }
}

static consolidated_subchunks_t *create_consolidated_subchunk(const int ndims)
{
    consolidated_subchunks_t *c=new consolidated_subchunks_t;
    c->start = (nc_size_t *)malloc(ndims*sizeof(nc_size_t));
    c->count = (nc_size_t *)malloc(ndims*sizeof(nc_size_t));
    c->ndims = ndims;

    return(c);
}

static consolidated_subchunks_t *consolidate_subchunk_iterval(const int        first,
                                                              const int        last,   /* inclusive */
                                                              const start_offset_length_service_t *sols,
                                                              const int        ndims,
                                                              const size_t    *dimlens,
                                                              const int        highest_contig_dim,
                                                              const nc_size_t *superchunk_start,
                                                              const nc_size_t *superchunk_count)
{
    if (first > last) {
        log_debug(netcdf_debug_level,"ERROR - first > last");
        return(NULL);
    }
    if (first == last) {
        consolidated_subchunks_t *c=create_consolidated_subchunk(ndims);
        memcpy(c->start, sols[first].start, ndims*sizeof(nc_size_t));
        memcpy(c->count, sols[first].count, ndims*sizeof(nc_size_t));
        nc_size_t dim_product=1;
        for (int i=0;i<ndims;i++) {
            dim_product *= (sols[first].start[i]-superchunk_start[i]);
        }
        log_debug(netcdf_debug_level,"dim_product(%ld)", dim_product);

        c->ndims=sols[first].ndims;
        c->nbytes=sols[first].length;

        c->offset_into_superchunk=(dim_product*sols[first].datatype_size);

        return(c);
    }
    for (int i=highest_contig_dim;i<ndims;i++) {
        if (sols[first].start[i] != sols[last].start[i]) {
            log_debug(netcdf_debug_level,"ERROR - discontiguous dimension don't have equal start coordinates "
                   "(sols[%d].start[%d](%ld) != sols[%d].start[%d](%ld)",
                   first, i, sols[first].start[i], last, i, sols[last].start[i]);
        }
    }
    int num_incomplete_dimensions=0;
    int *incomplete_dimensions=(int *)calloc(highest_contig_dim, sizeof(int));
    for (int i=0;i<highest_contig_dim;i++) {
        incomplete_dimensions[i]=-1;
        if ((sols[last].start[i]-sols[first].start[i]) != superchunk_count[i]-1) {
//            log_debug(netcdf_debug_level,"contiguous dimension doesn't span the entire superchunk "
//                   "sols[%d].start[%d](%ld)-(sols[%d].start[%d](%ld) != superchunk_count[%d]-1(%ld)",
//                   last, i, sols[last].start[i], first, i, sols[first].start[i], i, superchunk_count[i]-1);

            incomplete_dimensions[i]=(sols[last].start[i]-sols[first].start[i])+1;
            num_incomplete_dimensions++;
        }
    }
    consolidated_subchunks_t *c=NULL;
    if (num_incomplete_dimensions <= 1) {
        c=create_consolidated_subchunk(ndims);
        memcpy(c->start, sols[first].start, ndims*sizeof(nc_size_t));
        memcpy(c->count, superchunk_count, ndims*sizeof(nc_size_t));
        for (int i=0;i<highest_contig_dim;i++) {
            c->count[i] = (sols[last].start[i]-sols[first].start[i])+1;
            if (incomplete_dimensions[i] != -1) {
                c->count[i] = incomplete_dimensions[i];
            }
        }
        c->ndims=sols[first].ndims;
        c->nbytes=sols[first].length*(last-first+1);  // assume uniform length across all sols

        nc_size_t offset=0;
        for (int i=0;i<ndims;i++) {
            nc_size_t lower_dim_product=1;
            for (int j=i+1;j<ndims;j++) {
                lower_dim_product *= superchunk_count[j];
            }
//            log_debug(netcdf_debug_level, "lower_dim_product(%ld)", lower_dim_product);
            offset += ((sols[first].start[i]-superchunk_start[i]) * lower_dim_product);
        }
        c->offset_into_superchunk=(offset*sols[first].datatype_size);
    } else {
        log_debug(netcdf_debug_level,"there are multiple(%d) incomplete dimensions -- first subchunk is sols[%d].start[",
               num_incomplete_dimensions ,first);
        for (int k=0;k<ndims;k++) {
            log_debug(netcdf_debug_level,"%03ld,", sols[first].start[k]);
        }

//        log_debug(netcdf_debug_level,"]; last subchunk is sols[%d].start[", last);
//        for (int k=0;k<ndims;k++) {
//            log_debug(netcdf_debug_level,"%03ld,", sols[first].start[k]);
//        }
//        log_debug(netcdf_debug_level,"]");
    }
    free(incomplete_dimensions);

    return(c);
}

static consolidated_subchunks_map_t *consolidate_subchunks(const int        num_sols,
                                                           start_offset_length_service_t *sols,
                                                           const int        ndims,
                                                           const size_t    *dimlens,
                                                           const int        highest_contig_dim,
                                                           const nc_size_t *superchunk_start,
                                                           const nc_size_t *superchunk_count)
{
    /*
     * if sols was created by calc_offset_length, then sols is sorted by
     * sols.start[t][z][y][x] in ascending order.
     *
     * this also means sols.file_offset are in ascending order.  the rate
     * at which sols.file_offset increases is determined by:
     *   - (product of dims lower than highest_contig_dim)*sols.datatype_size
     *
     */

    consolidated_subchunks_map_t *map=new consolidated_subchunks_map_t;

    if (num_sols == 0) {
        // nothing at all.  return the empty map.

        return(map);
    }
    if (num_sols == 1) {
        // nothing to consolidate.

        consolidated_subchunks_vector_t *v=new consolidated_subchunks_vector_t;
        (*map)[sols[0].service]=v;
        consolidated_subchunks_t *c = consolidate_subchunk_iterval(0,
                                                                   0,
                                                                   sols,
                                                                   ndims,
                                                                   dimlens,
                                                                   highest_contig_dim,
                                                                   superchunk_start,
                                                                   superchunk_count);
        v->insert(v->end(), c);

        return(map);
    }



    // start with a sanity check
    nc_size_t lower_dim_product=1;
    for (int i=highest_contig_dim;i<ndims;i++) {
        lower_dim_product *= dimlens[i];
    }
    nc_size_t bytes_between_adjacent_subchunks=lower_dim_product*sols[0].datatype_size;
    if (bytes_between_adjacent_subchunks != (sols[1].offset-sols[0].offset)) {
        // way off
        log_debug(netcdf_debug_level,"file_offset calculation is not right.");
    } else {
        log_debug(netcdf_debug_level,"file_offset calculation is right on.");
    }




    int first_adjacent_chunk=0;
    for (int i=1;i<num_sols;i++) {
//        if ((sols[i].offset-sols[i-1].offset) == bytes_between_adjacent_subchunks) {
//            continue;
//        }
        if (sols[i].service == sols[i-1].service) {
            continue;
        }
        if (sols[first_adjacent_chunk].service != sols[i-1].service) {
            log_debug(netcdf_debug_level,"bummer these adjacent subchunks span multiple servers.  continue on for now.  FIX ME!!!!!!!!!!!!!!!!!");
        }
        consolidated_subchunks_t *c = consolidate_subchunk_iterval(first_adjacent_chunk,
                                                                   i-1,
                                                                   sols,
                                                                   ndims,
                                                                   dimlens,
                                                                   highest_contig_dim,
                                                                   superchunk_start,
                                                                   superchunk_count);

        consolidated_subchunks_vector_t *v=(*map)[sols[first_adjacent_chunk].service];
        if (v == NULL) {
            v=new consolidated_subchunks_vector_t;
            (*map)[sols[first_adjacent_chunk].service]=v;
        }
        v->insert(v->end(), c);

        first_adjacent_chunk=i;
    }
    consolidated_subchunks_t *c = consolidate_subchunk_iterval(first_adjacent_chunk,
            num_sols-1,
            sols,
            ndims,
            dimlens,
            highest_contig_dim,
            superchunk_start,
            superchunk_count);

    consolidated_subchunks_vector_t *v=(*map)[sols[first_adjacent_chunk].service];
    if (v == NULL) {
        v=new consolidated_subchunks_vector_t;
        (*map)[sols[first_adjacent_chunk].service]=v;
    }
    v->insert(v->end(), c);

    return(map);
}

consolidated_subchunks_map_t *netcdf_create_subchunks(const superchunk_t *chunk,
                                                      const int          *dimids,
                                                      const size_t       *dimlens,
                                                      const nc_size_t     bytes_per_server2)
{
    /*
     * I messed around alot before I came the realization that
     * unless an entire chunkwise Z plane fits on a server, you
     * have to spread subchunks of your chunk across the servers.
     * Guess what?  That's how mpich ADIO does it.
     *
     * This method ignores striding.  Fix it?
     */

    int rc=NC_NOERR;

    int unlimdimid=-1;
    int unlimdim_index=-1;

    start_offset_length_service_t *sols=NULL;
    int ol_index=0;
    nc_size_t *current_start=NULL;
    int highest_contig_dimension=0;
    int num_contig_subchunks=0;
    consolidated_subchunks_map_t *map=NULL;
    consolidated_subchunks_map_iterator_t map_iter;

    nc_size_t bytes_per_server;
    nc_size_t dim_product=1;
    log_debug(netcdf_debug_level, "ndims(%d)", chunk->ndims);
    for (int i=0;i<chunk->ndims;i++) {
        log_debug(netcdf_debug_level, "dimlens[%d](%d)", i, dimlens[i]);
        dim_product *= dimlens[i];
    }
    dim_product *= chunk->datatype_size;
    unlimdimid=-1;
    rc = nc_inq_unlimdim(chunk->ncid, &unlimdimid); /* get ID of unlimited dimension */
    if (rc != NC_NOERR) {
        log_error(netcdf_debug_level, "could not get unlimdimid");
        goto cleanup;
    }
    unlimdim_index=-1;
    if (unlimdimid != -1) {
        for (int i=0;i<chunk->ndims;i++) {
            if (unlimdimid == dimids[i]) {
                unlimdim_index=i;
                break;
            }
        }
        if ((unlimdim_index != 0) || (chunk->count[unlimdim_index] != 1)) {
            log_warn(netcdf_debug_level, "the unlimited dimension(%d) is not dim 0 or "
                     "the count in the unlimited dim(%ld) is not 1.  I don't know what to do.",
                     unlimdim_index, chunk->count[unlimdim_index]);
        } else if (chunk->start[unlimdim_index] > 0) {
            nc_size_t lower_dim_product=1;
            for (int i=1;i<chunk->ndims;i++) {
                lower_dim_product *= dimlens[i];
            }
            nc_size_t start_offset = (chunk->start[unlimdim_index]-1) * lower_dim_product * chunk->datatype_size;
            dim_product -= start_offset;
        }
    }
    bytes_per_server = dim_product/nc_cfg.num_servers;
    if (dim_product%nc_cfg.num_servers > 0) {
        bytes_per_server++;
    }
    log_debug(netcdf_debug_level, "dim_product(%lu) datatype_size(%d) num_servers(%d) dim_product/num_servers(%lu) bytes_per_server(%lu)",
            dim_product, chunk->datatype_size, nc_cfg.num_servers, dim_product/nc_cfg.num_servers, bytes_per_server);

    /* First, calculate the number of contiguous subchunks there are. */
    highest_contig_dimension=0;
    for (int i=chunk->ndims-1;i>=0;i--) {
        if (chunk->count[i] < dimlens[i]) {
            // the previous dim was the last contig dim
            highest_contig_dimension = i;
            break;
        }
    }
    num_contig_subchunks=1;
    for (int i=0;i<highest_contig_dimension;i++) {
        // calc the product of the discontig dims.  that is the number of
        // contiguous subchunks.
        num_contig_subchunks *= chunk->count[i];
    }
    /* alloc memory for the list of offset and lengths of the contiguous subchunks */
    sols=(start_offset_length_service_t *)calloc(num_contig_subchunks, sizeof(start_offset_length_service_t));
    ol_index=0;

    log_debug(netcdf_debug_level,"ndims(%ld) highest_contig_dimension(%d) num_contig_subchunks(%d)", chunk->ndims, highest_contig_dimension, num_contig_subchunks);

    double callTime;
    Start_Timer(callTime);
    current_start=(nc_size_t *)calloc(chunk->ndims, sizeof(nc_size_t));
    calc_offset_length(chunk->ndims,
                       dimlens,
                       highest_contig_dimension,
                       0,
                       0,
                       current_start,
                       sols,
                       &ol_index,
                       chunk->start,
                       chunk->count,
                       chunk->datatype_size);
    Stop_Timer("1st calc_offset_length", callTime);
//    for (int i=0;i<num_contig_subchunks;i++) {
//        log_debug(LOG_ALL,
//                "ncid(%d) varid(%d) start[%04ld,%04ld,%04ld,%04ld] count[%04ld,%04ld,%04ld,%04ld] file_offset(%07ld) length(%03ld)",
//                chunk->ncid, chunk->varid,
//                sols[i].start[0], sols[i].start[1], sols[i].start[2], sols[i].start[3],
//                sols[i].count[0], sols[i].count[1], sols[i].count[2], sols[i].count[3],
//                sols[i].offset, sols[i].length);
//    }

    Start_Timer(callTime);
    assign_service(dim_product, bytes_per_server, sols, num_contig_subchunks);
    Stop_Timer("assign_service", callTime);

    Start_Timer(callTime);
    map=consolidate_subchunks(num_contig_subchunks,
                              sols,
                              chunk->ndims,
                              dimlens,
                              highest_contig_dimension,
                              chunk->start,
                              chunk->count);
    Stop_Timer("consolidate_subchunks", callTime);
    Start_Timer(callTime);
    for (int i=0;i<num_contig_subchunks;i++) {
        log_debug(netcdf_debug_level,
                "ncid(%d) varid(%d) start[%04ld,%04ld,%04ld,%04ld] count[%04ld,%04ld,%04ld,%04ld] file_offset(%07ld) length(%03ld) service(%02d)",
                chunk->ncid, chunk->varid,
                sols[i].start[0], sols[i].start[1], sols[i].start[2], sols[i].start[3],
                sols[i].count[0], sols[i].count[1], sols[i].count[2], sols[i].count[3],
                sols[i].offset, sols[i].length, sols[i].service);
    }
    Stop_Timer("print consolidated subchunks", callTime);

    Start_Timer(callTime);
    map_iter = map->begin();
    for (;map_iter != map->end(); map_iter++) {
        consolidated_subchunks_vector_iterator_t vector_iter = (*map_iter).second->begin();
        for (;vector_iter != (*map_iter).second->end(); vector_iter++) {
            log_debug(netcdf_debug_level,
                    "ncid(%d) varid(%d) service(%02ld) super_start[%04ld,%04ld,%04ld,%04ld] start[%04ld,%04ld,%04ld,%04ld] count[%04ld,%04ld,%04ld,%04ld] offset_into_superchunk(%ld)",
                    chunk->ncid, chunk->varid,
                    (*map_iter).first,
                    chunk->start[0], chunk->start[1], chunk->start[2], chunk->start[3],
                    (*vector_iter)->start[0], (*vector_iter)->start[1], (*vector_iter)->start[2], (*vector_iter)->start[3],
                    (*vector_iter)->count[0], (*vector_iter)->count[1], (*vector_iter)->count[2], (*vector_iter)->count[3],
                    (*vector_iter)->offset_into_superchunk);
        }
    }
    Stop_Timer("print final subchunks", callTime);

cleanup:
    if (sols)          free(sols);
    if (current_start) free(current_start);

    return(map);
}
