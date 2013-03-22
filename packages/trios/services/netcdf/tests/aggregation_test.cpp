#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pnetcdf.h>
#include <mpi.h>

#include "Trios_nssi_client.h"
#include "netcdf_debug.h"

#include "aggregation.h"

int  ncid=0;                 /* netCDF id */
NNTI_peer_t caller;

int varid=1;
int ndims=4;

int dims[4] = {3, 3, 3, 3};
MPI_Offset *new_chunk_start=NULL;
MPI_Offset *new_chunk_count=NULL;
MPI_Offset *new_chunk_stride=NULL;

void add_chunk(int *buf)
{
    int buf_len;
    int num_elements;
    aggregation_chunk_details_t *chunk;

    new_chunk_count=(MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));
    new_chunk_stride=(MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));

    new_chunk_count[0]=dims[0];new_chunk_count[1]=dims[1];new_chunk_count[2]=dims[2];new_chunk_count[3]=dims[3];
    new_chunk_stride[0]=1;new_chunk_stride[1]=1;new_chunk_stride[2]=1;new_chunk_stride[3]=1;
    buf_len=dims[0]*dims[1]*dims[2]*dims[3]*sizeof(int);
    num_elements=dims[0]*dims[1]*dims[2]*dims[3];


    add_participant_for_file(ncid,
                             &caller,
                             WRITE_AGGREGATE_INDEPENDENT);

    chunk = new aggregation_chunk_details_t;
    chunk->ncid = ncid;
    chunk->varid = varid;
    chunk->ndims = ndims;

    chunk->buf = buf;

    chunk->atype = NC_INT;
    chunk->len   = buf_len;

    chunk->datatype      = MPI_INT;
    chunk->num_elements  = num_elements;

    chunk->start  = new_chunk_start;
    chunk->count  = new_chunk_count;
    chunk->stride = new_chunk_stride;

    add_participant_chunk(&caller, chunk);

}

/*
 * NB:  The dims of the individual chunks is the same as the dimensions
 * of the overall array in terms of chunks.  So if the dimensions of a chunk
 * is 3x3x3, then the overall array will have 3x3x3 chunks or 9x9x9 elements.
 */
int
main(int argc, char **argv)
{
    int test_chunk_count;
    int **chunk;

    int i=0,j=0,k=0,l=0;
    int buf_len;

    MPI_Init(&argc, &argv);

    netcdf_debug_level=LOG_ALL;
    /* initialize and enable logging */
    logger_init(netcdf_debug_level, NULL);

    test_chunk_count=1;
    buf_len=1;
    for (i=0;i<ndims;i++) {
        test_chunk_count *= dims[i];
        buf_len *= dims[i];
    }
    buf_len *= sizeof(int);

    chunk=(int **)malloc(test_chunk_count*sizeof(int *));
    for (i=0;i<test_chunk_count;i++) {
        chunk[i]=(int *)malloc(buf_len*sizeof(int));
        memset(chunk[i], 65+i, buf_len);
    }

    int chunk_index=0;
    for (i=0;i<dims[0];i++) {
        for (j=0;j<dims[1];j++) {
            for (k=0;k<dims[2];k++) {
                for (l=0;l<dims[2];l++) {
                    new_chunk_start=(MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));
                    new_chunk_start[0]=i*dims[0]; new_chunk_start[1]=j*dims[1]; new_chunk_start[2]=k*dims[2]; new_chunk_start[3]=l*dims[3];
                    add_chunk(chunk[chunk_index++]);
                }
            }
        }
    }

    while (try_aggregation(ncid, varid) == TRUE);

    int chunk_count;
    aggregation_chunk_details_t **chunks = get_chunks(ncid, varid, &chunk_count);
    for (i=0;i<chunk_count;i++) {
        print_chunk(chunks[i]);
    }
    free(chunks);

    cleanup_aggregation_chunks(ncid, varid);

    free(chunk);

    MPI_Finalize();

    return 0;
}
