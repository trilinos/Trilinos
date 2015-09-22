#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <map>
#include <list>
#include <vector>

using namespace std;

#define NDIMS 4

#define COUNT_DIM1 1
#define COUNT_DIM2 8
#define COUNT_DIM3 16
#define COUNT_DIM4 16

#define DATATYPE_SIZE 4

#define PROC_DIM1 1
#define PROC_DIM2 4
#define PROC_DIM3 4
#define PROC_DIM4 4

#define NUM_PROCS (PROC_DIM1*PROC_DIM2*PROC_DIM3*PROC_DIM4)
#define CLIENTS_PER_SERVER 4
#define BYTES_PER_SERVER (COUNT_DIM1*COUNT_DIM2*COUNT_DIM3*COUNT_DIM4*DATATYPE_SIZE*CLIENTS_PER_SERVER)

#define DIM1  (COUNT_DIM1*PROC_DIM1)
#define DIM2  (COUNT_DIM2*PROC_DIM2)
#define DIM3  (COUNT_DIM3*PROC_DIM3)
#define DIM4  (COUNT_DIM4*PROC_DIM4)

struct start_offset_length_service {
    long start[NDIMS];  /* coordinates of a contiguous subchunk */
    long count[NDIMS];  /* dimensions of a contiguous subchunk */
    long ndims;
    long offset;
    long length;
    long datatype_size;
    long service;
};
typedef struct start_offset_length_service start_offset_length_service_t;

struct consolidated_subchunks {
    long start[NDIMS];  /* coordinates of a consolidated subchunk */
    long count[NDIMS];  /* dimensions of a consolidated subchunk */
    long ndims;
    long offset_into_superchunk;  /* offset into the buffer that is the chunk that contains this subchunk */
};
typedef struct consolidated_subchunks consolidated_subchunks_t;

typedef list<consolidated_subchunks_t *> consolidated_subchunks_list_t;
typedef vector<consolidated_subchunks_t *> consolidated_subchunks_vector_t;
typedef vector<consolidated_subchunks_t *>::iterator consolidated_subchunks_vector_iterator_t;

typedef map<long, consolidated_subchunks_vector_t *> consolidated_subchunks_map_t;
typedef map<long, consolidated_subchunks_vector_t *>::iterator consolidated_subchunks_map_iterator_t;


void calc_offset_length(long ndims,
                        long *dimlens,
                        long highest_contig_dim,
                        long current_dim,
                        long file_offset,
                        long *current_index,
                        start_offset_length_service_t *sols,
                        long *ol_index,
                        long *start,
                        long *count,
                        long datatype_size)
{
    if (current_dim < highest_contig_dim) {
        current_index[current_dim]=start[current_dim];
        for (int i=0;i<count[current_dim];i++) {
            long my_offset_adder = current_index[current_dim];
            for (int j=current_dim+1;j<ndims;j++) {
                my_offset_adder *= dimlens[j];
            }
            my_offset_adder *= datatype_size;

            current_index[current_dim+1]=0;
            calc_offset_length(ndims, dimlens, highest_contig_dim, current_dim+1, file_offset+my_offset_adder,
                               current_index, sols, ol_index, start, count, datatype_size);
            current_index[current_dim]++;
        }
    } else {
        current_index[current_dim]=start[current_dim];
        long my_offset_adder = start[current_dim];
        for (int j=current_dim+1;j<ndims;j++) {
            my_offset_adder *= dimlens[j];
        }
        my_offset_adder *= datatype_size;
        long length=1;
        for (int j=current_dim;j<ndims;j++) {
            length *= count[j];
        }
        length *= datatype_size;

        // populate the sols for this subchunk
        memcpy(sols[*ol_index].start, current_index, sizeof(sols[*ol_index].start));
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


//        printf("start[");
//        for (int k=0;k<ndims;k++) {
//            printf("%03lu,", current_index[k]);
//        }
//        printf("] file_offset(%lu) length(%lu)\n", file_offset+my_offset_adder, length);
    }
}


void assign_service(long bytes_per_service, start_offset_length_service_t *sols, long num_sols)
{
    for(long i=0;i<num_sols;i++) {
        sols[i].service = sols[i].offset/bytes_per_service;
    }
}

consolidated_subchunks_t *consolidate_subchunk_iterval(long first,
                                                       long last,   /* inclusive */
                                                       start_offset_length_service_t *sols,
                                                       long ndims,
                                                       long *dimlens,
                                                       long highest_contig_dim,
                                                       long *superchunk_start,
                                                       long *superchunk_count)
{
    if (first > last) {
        printf("ERROR - first > last");
        return(NULL);
    }
    if (first == last) {
        consolidated_subchunks_t *c=new consolidated_subchunks_t;
        memcpy(c->start, sols[first].start, sols[first].ndims*sizeof(long));
        memcpy(c->count, sols[first].count, sols[first].ndims*sizeof(long));
        c->ndims=sols[first].ndims;
        long dim_product=1;
        for (int i=0;i<ndims;i++) {
            dim_product *= (sols[first].start[i]-superchunk_start[i]);
        }
        printf("dim_product(%ld)\n", dim_product);
        c->offset_into_superchunk=(dim_product*sols[first].datatype_size);

        return(c);
    }
    for (long i=highest_contig_dim;i<ndims;i++) {
        if (sols[first].start[i] != sols[last].start[i]) {
            printf("ERROR - discontiguous dimension don't have equal start coordinates "
                   "(sols[%ld].start[%ld](%ld) != sols[%ld].start[%ld](%ld)\n",
                   first, i, sols[first].start[i], last, i, sols[last].start[i]);
        }
    }
    long num_incomplete_dimensions=0;
    long *incomplete_dimensions=(long *)calloc(highest_contig_dim, sizeof(long));
    for (long i=0;i<highest_contig_dim;i++) {
        incomplete_dimensions[i]=-1;
        if ((sols[last].start[i]-sols[first].start[i]) != superchunk_count[i]-1) {
            printf("contiguous dimension doesn't span the entire superchunk "
                   "sols[%ld].start[%ld](%ld)-(sols[%ld].start[%ld](%ld) != superchunk_count[%ld]-1(%ld)\n",
                   last, i, sols[last].start[i], first, i, sols[first].start[i], i, superchunk_count[i]-1);

            incomplete_dimensions[i]=(sols[last].start[i]-sols[first].start[i])+1;
            num_incomplete_dimensions++;
        }
    }
    consolidated_subchunks_t *c=NULL;
    if (num_incomplete_dimensions <= 1) {
        c=new consolidated_subchunks_t;
        memcpy(c->start, sols[first].start, sols[first].ndims*sizeof(long));
        memcpy(c->count, superchunk_count, sols[first].ndims*sizeof(long));
        for (long i=0;i<highest_contig_dim;i++) {
            c->count[i] = (sols[last].start[i]-sols[first].start[i])+1;
            if (incomplete_dimensions[i] != -1) {
                c->count[i] = incomplete_dimensions[i];
            }
        }
        c->ndims=sols[first].ndims;
        long offset=0;
        for (int i=0;i<ndims;i++) {
            long lower_dim_product=1;
            for (int j=i+1;j<ndims;j++) {
//                lower_dim_product *= ((sols[first].start[j]-superchunk_start[j]+1) * superchunk_count[j]);
                lower_dim_product *= superchunk_count[j];
            }
            printf("lower_dim_product(%ld)\n", lower_dim_product);
            offset += ((sols[first].start[i]-superchunk_start[i]) * lower_dim_product);
        }
        c->offset_into_superchunk=(offset*sols[first].datatype_size);
    } else {
        printf("there are multiple(%ld) incomplete dimensions -- first subchunk is sols[%ld].start[",
               num_incomplete_dimensions ,first);
        for (int k=0;k<ndims;k++) {
            printf("%03lu,", sols[first].start[k]);
        }

        printf("]; last subchunk is sols[%ld].start[", last);
        for (int k=0;k<ndims;k++) {
            printf("%03lu,", sols[first].start[k]);
        }
        printf("]\n");
    }

//    for (long i=first;i<=last;i++) {
//
//    }

    return(c);
}

consolidated_subchunks_map_t *consolidate_subchunks(long num_sols,
                                                    start_offset_length_service_t *sols,
                                                    long ndims,
                                                    long *dimlens,
                                                    long highest_contig_dim,
                                                    long *superchunk_start,
                                                    long *superchunk_count)
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
    long lower_dim_product=1;
    for (int i=highest_contig_dim;i<ndims;i++) {
        lower_dim_product *= dimlens[i];
    }
    long bytes_between_adjacent_subchunks=lower_dim_product*sols[0].datatype_size;
    if (bytes_between_adjacent_subchunks != (sols[1].offset-sols[0].offset)) {
        // way off
        printf("file_offset calculation is not right.\n");
    } else {
        printf("file_offset calculation is right on.\n");
    }




    long first_adjacent_chunk=0;
    for (long i=0;i<num_sols;i++) {
        if ((sols[i+1].offset-sols[i].offset) == bytes_between_adjacent_subchunks) {
            continue;
        }
        if (sols[first_adjacent_chunk].service != sols[i].service) {
            printf("bummer these adjacent subchunks span multiple servers.  continue on for now.  FIX ME!!!!!!!!!!!!!!!!!\n");
        }
        consolidated_subchunks_t *c = consolidate_subchunk_iterval(first_adjacent_chunk,
                                                                   i,
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

        first_adjacent_chunk=i+1;
    }

    return(map);
}

int main(int argc, char **argv)
{
        long ndims=NDIMS;
        long dimlens[NDIMS]={DIM1, DIM2, DIM3, DIM4};
        long start[NDIMS];
        long count[NDIMS]={COUNT_DIM1,COUNT_DIM2,COUNT_DIM3,COUNT_DIM4};
        long num_servers;

        num_servers = (NUM_PROCS/CLIENTS_PER_SERVER);
        if (NUM_PROCS%CLIENTS_PER_SERVER > 0) num_servers++;

        printf("%lu servers are required to service %d clients\n", num_servers, NUM_PROCS);

        for (long a=0;a<DIM1;a+=COUNT_DIM1) {
            for (long b=0;b<DIM2;b+=COUNT_DIM2) {
                for (long c=0;c<DIM3;c+=COUNT_DIM3) {
                    for (long d=0;d<DIM4;d+=COUNT_DIM4) {
                        start[0]=a;
                        start[1]=b;
                        start[2]=c;
                        start[3]=d;


/* ************************************************************************* */
                        /*
                         * I messed around alot before I came the realization that
                         * unless an entire chunkwise Z plane fits on a server, you
                         * have to spread subchunks of your chunk across the servers.
                         * Guess what?  That's how mpich ADIO does it.
                         *
                         * This method ignores striding.  Fix it?
                         */

                        /* First, calculate the number of contiguous subchunks there are. */
                        long highest_contig_dimension=0;
                        for (int i=ndims-1;i>=0;i--) {
                            if (count[i] < dimlens[i]) {
                                // the previous dim was the last contig dim
                                highest_contig_dimension = i;
                                break;
                            }
                        }
                        long num_contig_subchunks=1;
                        for (int i=0;i<highest_contig_dimension;i++) {
                            // calc the product of the discontig dims.  that is the number of
                            // contiguous subchunks.
                            num_contig_subchunks *= count[i];
                        }
                        /* alloc memory for the list of offset and lengths of the contiguous subchunks */
                        start_offset_length_service_t *sols=
                            (start_offset_length_service_t *)calloc(num_contig_subchunks, sizeof(start_offset_length_service_t));
                        long ol_index=0;

                        printf("highest_contig_dimension(%ld) num_contig_subchunks(%ld)\n", highest_contig_dimension, num_contig_subchunks);

                        long *current_index=(long *)calloc(ndims, sizeof(long));
                        calc_offset_length(ndims,
                                           dimlens,
                                           highest_contig_dimension,
                                           0,
                                           0,
                                           current_index,
                                           sols,
                                           &ol_index,
                                           start,
                                           count,
                                           DATATYPE_SIZE);

                        assign_service(BYTES_PER_SERVER, sols, num_contig_subchunks);

                        consolidated_subchunks_map_t *map=consolidate_subchunks(num_contig_subchunks,
                                                                                sols,
                                                                                ndims,
                                                                                dimlens,
                                                                                highest_contig_dimension,
                                                                                start,
                                                                                count);

                        for (long i=0;i<num_contig_subchunks;i++) {
                            printf("start[%04ld,%04ld,%04ld,%04ld] count[%04ld,%04ld,%04ld,%04ld] file_offset(%07ld) length(%03ld) service(%02ld)\n",
                                    sols[i].start[0], sols[i].start[1], sols[i].start[2], sols[i].start[3],
                                    sols[i].count[0], sols[i].count[1], sols[i].count[2], sols[i].count[3],
                                    sols[i].offset, sols[i].length, sols[i].service);
                        }

                        consolidated_subchunks_map_iterator_t map_iter = map->begin();
                        for (;map_iter != map->end(); map_iter++) {
                            consolidated_subchunks_vector_iterator_t vector_iter = (*map_iter).second->begin();
                            for (;vector_iter != (*map_iter).second->end(); vector_iter++) {
                                printf("service(%02ld) super_start[%04ld,%04ld,%04ld,%04ld] start[%04ld,%04ld,%04ld,%04ld] count[%04ld,%04ld,%04ld,%04ld] offset_into_superchunk(%ld)\n",
                                        (*map_iter).first,
                                        start[0], start[1], start[2], start[3],
                                        (*vector_iter)->start[0], (*vector_iter)->start[1], (*vector_iter)->start[2], (*vector_iter)->start[3],
                                        (*vector_iter)->count[0], (*vector_iter)->count[1], (*vector_iter)->count[2], (*vector_iter)->count[3],
                                        (*vector_iter)->offset_into_superchunk);

                            }

                        }




/* ************************************************************************* */



/* ************************************************************************* */
//                        /*
//                         * (z+1)*y/CLIENT_PER_SERVER
//                         */
//                        unsigned long dim_product=1;
//                        unsigned long count_product=1;
//                        unsigned long servers_required=0;
//                        if ((ndims == 4) && (count[0]) == 1) {
//                            // assume dim1 is the timestep dimension
//                            dim_product=dimlens[1]*dimlens[2];
//                            count_product=count[1]*count[2];
//                            servers_required=(dim_product/count_product)/CLIENTS_PER_SERVER;
//                            if (((dim_product/count_product)%CLIENTS_PER_SERVER) > 0) {
//                                servers_required++;
//                            }
//                            if (servers_required > num_servers) {
//                                printf("(dim_product(%lu)/count_product(%lu)/CLIENTS_PER_SERVER(%ul))(%lu) > num_servers(%lu)  "
//                                        "you need more servers.\n",
//                                        dim_product, count_product, (dim_product/count_product/CLIENTS_PER_SERVER), num_servers);
//                            }
//                            svc_index =
//                        }
/* ************************************************************************* */

/* ************************************************************************* */
                        /* this method of server selection is for a very specific array/client layout.
                         * assumptions:
                         *  - ndims is 3 or 4.
                         *  - if ndims == 4, then
                         *    - dim1 is the timestep dimension.  each client only works on one timestep
                         *      at a time.  therefore, dim1 may be >= 0, but count[0] is always 1.
                         *    - dim2 is the Z dimension.  we can split here.
                         *    - dim3 is the Y dimension.  we can split here.
                         *    - dim4 is the X dimension.  we cannot split processors on this dimension
                         *      and get any decent aggregation.
                         *  - if ndims == 3, then
                         *    - if count[0] == 1, then
                         *      - dim1 is the timestep dimension.  each client only works on one timestep
                         *        at a time.  therefore, dim1 may be >= 0, but count[0] is always 1.
                         *      - dim2 is the Y dimension.  we can split here.
                         *      - dim3 is the X dimension.  we cannot split processors on this dimension
                         *        and get any decent aggregation.
                         *    - if count[0] > 1, then
                         *      - assume dim1 is the Z dimension.  we can split here.
                         *      - dim2 is the Y dimension.  we can split here.
                         *      - dim3 is the X dimension.  we cannot split processors on this dimension
                         *        and get any decent aggregation.
                         */

//                        unsigned long dim_product=1;
//                        unsigned long count_product=1;
//                        long split_dim=-1;
//                        if ((ndims == 4) && (count[0]) == 1) {
//                            dim_product=dimlens[1];
//                            count_product=count[1];
//                            if ((dim_product/count_product) >= num_servers) {
//                                split_dim=1;
//                            } else {
//                                dim_product *= dimlens[2];
//                                count_product *= count[2];
//                                if ((dim_product/count_product) >= num_servers) {
//                                    split_dim=2;
//                                }
//                            }
//                        }
/* ************************************************************************* */


/* ************************************************************************* */
//                        unsigned long dim_product=1;
//                        unsigned long count_product=1;
//                        long split_dim=-1;
//                        for (unsigned long i=0;i<ndims;i++) {
//                            dim_product *= dimlens[i];
//                            count_product *= count[i];
//                            if (dim_product/count_product >= num_servers) {
//                                // we have enough dimensions to divide up the data
//                                split_dim=i;
//                                break;
//                            }
//                        }
//                        if (split_dim == -1) {
//                            if ((a%COUNT_DIM == 0) &&
//                                (b%COUNT_DIM == 0) &&
//                                (c%COUNT_DIM == 0) &&
//                                (d%COUNT_DIM == 0)) {
//                                printf("(dim_product(%lu)/count_product(%lu))(%lu) < num_servers(%lu)  "
//                                        "using default_server is your best bet.\n",
//                                        dim_product, count_product, (dim_product/count_product), num_servers);
//                            }
//                        }
//                        unsigned long lower_dim_product=1;
//                        unsigned long start_offset=0;
//                        for (unsigned long i=0;i<split_dim+1;i++) {
//                            lower_dim_product=1;
//                            for (unsigned long j=i+1;j<split_dim;j++) {
//                                lower_dim_product *= dimlens[j];
//                            }
//                            start_offset += start[i]*lower_dim_product;
////                            printf("start_offset==%lu\n", start_offset);
//                        }
//                        unsigned long bytes_per_server=dim_product/num_servers;
//
//                        unsigned long svc_index=(start_offset/bytes_per_server);
//                        if ((a%COUNT_DIM == 0) &&
//                            (b%COUNT_DIM == 0) &&
//                            (c%COUNT_DIM == 0) &&
//                            (d%COUNT_DIM == 0)) {
//                            printf("start[%03ld,%03ld,%03ld,%03ld] split_dim(%lu) dim_product(%lu) start_offset(%lu) bytes_per_server(%lu) count_product(%lu) svc_index(%lu)\n",
//                                    a, b, c, d, split_dim, dim_product, start_offset, bytes_per_server, count_product, svc_index);
//                        }
/* ************************************************************************* */

/* ************************************************************************* */
//                        unsigned long dim_product=1;
//                        unsigned long bytes_per_server=0;
//                        unsigned long my_file_offset=0;
//                        unsigned long lower_dim_product=1;
//                        for (unsigned long i=0;i<ndims;i++) {
//                            dim_product *= dimlens[i];
//
//                            lower_dim_product=1;
//                            for (unsigned long j=i+1;j<ndims;j++) {
//                                lower_dim_product *= dimlens[j];
//                            }
//                            my_file_offset += start[i]*lower_dim_product;
//                            printf("lower_dim_product(%lu)\n", lower_dim_product);
//                            printf("start[%lu](%lu)\n", i, start[i]);
//                            printf("my_file_offset(%lu)\n", my_file_offset);
//                        }
//                        printf("=======================================\n");
//                        bytes_per_server=dim_product/num_servers;
//                        unsigned long svc_index=(my_file_offset/bytes_per_server);
//                        printf("start[%03ld,%03ld,%03ld,%03ld] dim_product(%lu) bytes_per_server(%lu) my_file_offset(%lu) svc_index(%lu)\n",
//                                a, b, c, d, dim_product, bytes_per_server, my_file_offset, svc_index);

/* ************************************************************************* */


/* ************************************************************************* */

//                        unsigned long svc_index;
//                        unsigned long dim_product=1;
//                        unsigned long elements_per_server=0;
//
//                        for (unsigned long i=0;i<ndims;i++) {
//                            dim_product *= dimlens[i];
//                            elements_per_server = dim_product/num_servers;
//                            if (elements_per_server > 0) {
//                                if (dim_product%num_servers > 0) {
//                                    elements_per_server++;
//                                }
//                                break;
//                            }
//                        }
//                        if (elements_per_server == 0) {
//                            printf("elements_per_server==0.  use default_svc.\n");
//                        } else {
//                            unsigned long my_file_offset=0;
//                            unsigned long lower_dim_product=1;
//                            for (unsigned long i=0;i<ndims;i++) {
//                                lower_dim_product=1;
//                                for (unsigned long j=i+1;j<ndims;j++) {
//                                    lower_dim_product *= dimlens[j];
//                                }
//                                my_file_offset += start[i]*lower_dim_product;
//                        if ((a%COUNT_DIM == 0) &&
//                            (b%COUNT_DIM == 0) &&
//                            (c%COUNT_DIM == 0) &&
//                            (d%COUNT_DIM == 0)) {
//                                    printf("lower_dim_product(%lu)\n", lower_dim_product);
//                                    printf("start[%lu](%lu)\n", i, start[i]);
//                                    printf("my_file_offset(%lu)\n", my_file_offset);
//                                }
//                            }
//                            svc_index = (my_file_offset/elements_per_server);
//
//                        if ((a%COUNT_DIM == 0) &&
//                            (b%COUNT_DIM == 0) &&
//                            (c%COUNT_DIM == 0) &&
//                            (d%COUNT_DIM == 0)) {
//                                printf("start[%03ld,%03ld,%03ld,%03ld] elements_per_server(%lu) my_file_offset(%lu) svc_index(%lu)\n",
//                                    a, b, c, d, elements_per_server, my_file_offset, svc_index);
//                            }
//                        }
/* ************************************************************************* */
                    }
                }
            }
        }
}
