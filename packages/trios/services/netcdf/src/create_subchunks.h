/*
 * create_subchunks.h
 *
 *  Created on: Mar 26, 2009
 *      Author: thkorde
 */

#ifndef CREATE_SUBCHUNKS_H_
#define CREATE_SUBCHUNKS_H_


#include <algorithm>
#include <map>
#include <list>
#include <vector>

using namespace std;

#include "Trios_nssi_types.h"
#include "netcdf_args.h"

#ifdef __cplusplus
extern "C" {
#endif

struct superchunk_t {
    int ncid;
    int varid;

    int ndims;

    void *buf; /* the data */

    size_t len;   /* length of buf in bytes */

    nc_type      buftype;  /* the netcdf type of the data in buf */
    size_t datatype_size;

    nc_size_t *start;  /* starting corner (eg. 0,0,0 is the origin of a cube) */
    nc_size_t *count;  /* num elements in each dimension (eg. 3,3,3 is a cube of size 3) */
    nc_size_t *stride; /* sampling interval (eg. 2,2,2 is every other element in each dimension) */
};
typedef struct superchunk_t superchunk_t;

struct consolidated_subchunks {
    nc_size_t *start;  /* coordinates of a consolidated subchunk */
    nc_size_t *count;  /* dimensions of a consolidated subchunk */
    int        ndims;
    int        nbytes;
    nc_size_t  offset_into_superchunk;  /* offset into the buffer that is the chunk that contains this subchunk */
};
typedef struct consolidated_subchunks consolidated_subchunks_t;

typedef list<consolidated_subchunks_t *> consolidated_subchunks_list_t;
typedef vector<consolidated_subchunks_t *> consolidated_subchunks_vector_t;
typedef vector<consolidated_subchunks_t *>::iterator consolidated_subchunks_vector_iterator_t;

typedef map<nc_size_t, consolidated_subchunks_vector_t *> consolidated_subchunks_map_t;
typedef map<nc_size_t, consolidated_subchunks_vector_t *>::iterator consolidated_subchunks_map_iterator_t;



consolidated_subchunks_map_t *netcdf_create_subchunks(const superchunk_t *chunk,
                                                      const int          *dimids,
                                                      const size_t       *dimlens,
                                                      const nc_size_t     bytes_per_server);

#ifdef __cplusplus
}
#endif

#endif /* CREATE_SUBCHUNKS_H_ */
