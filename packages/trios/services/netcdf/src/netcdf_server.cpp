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
/**  @file main.c
 *
 *   @brief Driver for the NSSI name server.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 *   $Revision: 1264 $.
 *   $Date: 2007-02-27 15:30:26 -0700 (Tue, 27 Feb 2007) $.
 */


#include "Trios_config.h"

#ifdef HAVE_TRIOS_PNETCDF

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif

#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <sys/mman.h>

#include <Trios_nssi_server.h>  /* for nssi service functions */
#include <Trios_logger.h>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_oblackholestream.hpp"

#include "pnetcdf.h"
#include <mpi.h>
#include <algorithm>

using namespace std;

#include "netcdf_args.h"
#include "netcdf_config.h"
#include "netcdf_debug.h"  /* netcdf_debug_level */

#include "aggregation.h"

#include "io_timer.h"


typedef char NNTI_url[NNTI_URL_LEN];


#ifdef __LIBCATAMOUNT__
#define ntohs(value) 0
#endif


static double create_time=0.0;
static double open_time  =0.0;
static double enddef_time=0.0;
static double close1_time=0.0;
static double close2_time=0.0;


/* -------------------- PRIVATE FUNCTIONS ---------- */


static int do_put_vars(int ncid, int varid,
               const MPI_Offset start[],
               const MPI_Offset count[],
               const MPI_Offset stride[],
               const void *buf, MPI_Offset bufcount,
               MPI_Datatype datatype)
{
    int rc=NC_NOERR;

    log_debug(netcdf_debug_level, "calling ncmpi_put_vars with direct data");
    rc = ncmpi_put_vars(ncid, varid, start, count, stride,
            buf, bufcount, datatype);

    return rc;
}

static int do_put_vars(aggregation_chunk_details_t **chunks, const int chunk_count)
{
    int rc=NC_NOERR;
    double call_time;

    for (int i=0;i<chunk_count;i++) {
        log_debug(netcdf_debug_level, "calling ncmpi_put_vars with agg data");
        Start_Timer(call_time);
        rc = ncmpi_put_vars(chunks[i]->ncid,
                            chunks[i]->varid,
                            chunks[i]->start,
                            chunks[i]->count,
                            chunks[i]->stride,
                            chunks[i]->buf,
                            chunks[i]->num_elements,
                            chunks[i]->datatype);
        Stop_Timer("ncmpi_put_vars", call_time);
        if (rc != NC_NOERR) {
            log_error(netcdf_debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
            goto cleanup;
        }
    }

cleanup:
    return rc;
}

static int do_put_vars_all(const int ncid, const int varid, aggregation_chunk_details_t **chunks, int chunk_count)
{
    int rc=NC_NOERR;
    int min_chunk_count=0;
    int max_chunk_count=0;

    /*
     *  max_chunk_count is the largest number of chunks held by any server.  this will
     *  determine is anyone has chunks for this varid.
     */
    MPI_Allreduce(&chunk_count, &max_chunk_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if (max_chunk_count == 0) {
        // no one has any chunks for this varid.  bail out.
        return NC_NOERR;
    }

    /*
     *  min_chunk_count is the smallest number of chunks held by any server.  this is
     *  the maximum number of collective calls allowed.  remaining chunks will be put
     *  in independent data mode.
     */
    MPI_Allreduce(&chunk_count, &min_chunk_count, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    for (int i=0;i<min_chunk_count;i++) {
        log_debug(netcdf_debug_level, "calling ncmpi_put_vars_all with agg data");

        log_debug(netcdf_debug_level,
                "ncid(%d) varid(%d) ndims(%d) start[%04ld,%04ld,%04ld,%04ld] count[%04ld,%04ld,%04ld,%04ld] num_elements(%d)",
                chunks[i]->ncid, chunks[i]->varid, chunks[i]->ndims,
                chunks[i]->start[0], chunks[i]->start[1], chunks[i]->start[2], chunks[i]->start[3],
                chunks[i]->count[0], chunks[i]->count[1], chunks[i]->count[2], chunks[i]->count[3],
                chunks[i]->num_elements);
        rc = ncmpi_put_vars_all(chunks[i]->ncid,
                chunks[i]->varid,
                chunks[i]->start,
                chunks[i]->count,
                chunks[i]->stride,
                chunks[i]->buf,
                chunks[i]->num_elements,
                chunks[i]->datatype);

        log_debug(netcdf_debug_level, "ncmpi_put_vars_all: rc==%d", rc);
        if (rc != NC_NOERR) {
            log_error(netcdf_debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
            goto cleanup;
        }
    }

    if (min_chunk_count < chunk_count) {
        ncmpi_begin_indep_data(ncid);
        for (int i=min_chunk_count;i<chunk_count;i++) {
            log_debug(netcdf_debug_level, "calling ncmpi_put_vars with agg data");
            rc = ncmpi_put_vars(chunks[i]->ncid,
                    chunks[i]->varid,
                    chunks[i]->start,
                    chunks[i]->count,
                    chunks[i]->stride,
                    chunks[i]->buf,
                    chunks[i]->num_elements,
                    chunks[i]->datatype);

            log_debug(netcdf_debug_level, "ncmpi_put_vars: rc==%d", rc);
            if (rc != NC_NOERR) {
                log_error(netcdf_debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
                goto cleanup;
            }
        }
        ncmpi_end_indep_data(ncid);
    }

cleanup:
    return rc;
}

static int flush_aggregated_chunks(
        const int ncid,
        const int varid)
{
    int rc=NC_NOERR;
    int chunk_count=0;
    double callTime;

    Start_Timer(callTime);
    try_aggregation(ncid, varid);
    Stop_Timer("agg", callTime);

    aggregation_chunk_details_t **chunks = get_chunks(ncid, varid, &chunk_count);
    if (chunk_count == 0) {
        log_debug(netcdf_debug_level, "chunk_count is 0, but I must participate in the collective.");
    } else {
        log_debug(netcdf_debug_level, "chunk_count is %d.", chunk_count);
    }
//    log_level old_log=netcdf_debug_level;
//    netcdf_debug_level = netcdf_debug_level;
//    for (int i=0;i<chunk_count;i++) {
//    	print_chunk(chunks[i]);
//    }
//    netcdf_debug_level = old_log;

    if (use_collective(ncid)) {
        Start_Timer(callTime);
        rc = do_put_vars_all(ncid, varid, chunks, chunk_count);
        if (rc == NC_EINDEP) {
            rc = do_put_vars(chunks, chunk_count);
        }
        Stop_Timer("put_vars_all", callTime);

        free(chunks);
        cleanup_aggregation_chunks(ncid, varid);

        log_debug(netcdf_debug_level, "do_put_vars_all: rc==%d", rc);
        if (rc != NC_NOERR) {
            log_error(netcdf_debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
        }
    } else {
        Start_Timer(callTime);
        rc = do_put_vars(chunks, chunk_count);
        if (rc == NC_ENOTINDEP) {
            rc = do_put_vars_all(ncid, varid, chunks, chunk_count);
        }
        Stop_Timer("put_vars", callTime);

        free(chunks);
        cleanup_aggregation_chunks(ncid, varid);

        log_debug(netcdf_debug_level, "do_put_vars: rc==%d", rc);
        if (rc != NC_NOERR) {
            log_error(netcdf_debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
        }
    }

    return rc;
}
static int flush_aggregated_chunks(const int ncid)
{
    int rc=NC_NOERR;
    int chunk_count=0;
    double callTime;

    if (use_aggregation(ncid)) {
        log_debug(netcdf_debug_level, "flushing aggregation chunks");

        if (use_collective(ncid)) {
            log_debug(netcdf_debug_level, "using collective");
            // get num vars
            // iter thru vars
            //  - try_agg on varid
            //  - get_chunks
            //  - do_put_vars_all

            int var_count=0;
            int varid=0;
            int chunk_count=0;

            /* get varids for this group */
            if ((rc = ncmpi_inq_nvars(ncid, &var_count)) != NC_NOERR) {
                log_error(netcdf_debug_level, "Could not get varids: %s", ncmpi_strerror(rc));
                goto cleanup;
            }
            for (int i=0;i<var_count;i++) {
                varid=i;

                flush_aggregated_chunks(ncid, varid);
            }

        } else {
            log_debug(netcdf_debug_level, "using independent");
            /* this is the easy one.  use independent mode. */

            Start_Timer(callTime);
            try_aggregation(ncid);  // aggregate all varids in this file
            Stop_Timer("agg", callTime);

            aggregation_chunk_details_t **chunks = get_chunks(ncid, &chunk_count);
            if (chunk_count == 0) {
                log_debug(netcdf_debug_level, "chunk_count is 0.  nothing to do.");
                goto cleanup;
            } else {
                log_debug(netcdf_debug_level, "chunk_count is %d.", chunk_count);
            }
//            log_level old_log=netcdf_debug_level;
//            netcdf_debug_level = netcdf_debug_level;
//            for (int i=0;i<chunk_count;i++) {
//                print_chunk(chunks[i]);
//            }
//            netcdf_debug_level = old_log;

            Start_Timer(callTime);
            rc = do_put_vars(chunks, chunk_count);
            Stop_Timer("put_vars", callTime);

            free(chunks);
            cleanup_aggregation_chunks(ncid);

            if (rc != NC_NOERR) {
                log_error(netcdf_debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
            }
        }
    }

cleanup:

    return rc;
}

static int flush_cached_chunks(
        const int ncid,
        const int varid)
{
    int rc=NC_NOERR;
    int chunk_count=0;
    double callTime;

    aggregation_chunk_details_t **chunks = get_chunks(ncid, varid, &chunk_count);
    if (chunk_count == 0) {
        log_debug(netcdf_debug_level, "chunk_count is 0, but I must participate in the collective.");
    } else {
        log_debug(netcdf_debug_level, "chunk_count is %d.", chunk_count);
    }
//    log_level old_log=netcdf_debug_level;
//    netcdf_debug_level = netcdf_debug_level;
//    for (int i=0;i<chunk_count;i++) {
//    	print_chunk(chunks[i]);
//    }
//    netcdf_debug_level = old_log;

    if (use_collective(ncid)) {
        Start_Timer(callTime);
        rc = do_put_vars_all(ncid, varid, chunks, chunk_count);
        if (rc == NC_EINDEP) {
            rc = do_put_vars(chunks, chunk_count);
        }
        Stop_Timer("put_vars_all", callTime);

        free(chunks);
        cleanup_aggregation_chunks(ncid, varid);

        if (rc != NC_NOERR) {
            log_error(netcdf_debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
        }
    } else {
        Start_Timer(callTime);
        rc = do_put_vars(chunks, chunk_count);
        if (rc == NC_ENOTINDEP) {
            rc = do_put_vars_all(ncid, varid, chunks, chunk_count);
        }
        Stop_Timer("put_vars", callTime);

        free(chunks);
        cleanup_aggregation_chunks(ncid, varid);

        if (rc != NC_NOERR) {
            log_error(netcdf_debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
        }
    }

    return rc;
}
static int flush_cached_chunks(const int ncid)
{
    int rc=NC_NOERR;
    int chunk_count=0;
    double callTime;

    if (use_caching(ncid)) {
        log_debug(netcdf_debug_level, "flushing cached chunks");

        if (use_collective(ncid)) {
            log_debug(netcdf_debug_level, "using collective");
            // get num vars
            // iter thru vars
            //  - try_agg on varid
            //  - get_chunks
            //  - do_put_vars_all

            int var_count=0;
            int varid=0;
            int chunk_count=0;

            /* get varids for this group */
            if ((rc = ncmpi_inq_nvars(ncid, &var_count)) != NC_NOERR) {
                log_error(netcdf_debug_level, "Could not get varids: %s", ncmpi_strerror(rc));
                goto cleanup;
            }
            for (int i=0;i<var_count;i++) {
                varid=i;

                flush_cached_chunks(ncid, varid);
            }

        } else {
            log_debug(netcdf_debug_level, "using independent");
            /* this is the easy one.  use independent mode. */

            aggregation_chunk_details_t **chunks = get_chunks(ncid, &chunk_count);
            if (chunk_count == 0) {
                log_debug(netcdf_debug_level, "chunk_count is 0.  nothing to do.");
                goto cleanup;
            } else {
                log_debug(netcdf_debug_level, "chunk_count is %d.", chunk_count);
            }
//            log_level old_log=netcdf_debug_level;
//            netcdf_debug_level = netcdf_debug_level;
//            for (int i=0;i<chunk_count;i++) {
//                print_chunk(chunks[i]);
//            }
//            netcdf_debug_level = old_log;

            Start_Timer(callTime);
            rc = do_put_vars(chunks, chunk_count);
            Stop_Timer("put_vars", callTime);

            free(chunks);
            cleanup_aggregation_chunks(ncid);

            if (rc != NC_NOERR) {
                log_error(netcdf_debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
                goto cleanup;
            }
        }
    }

cleanup:

    return rc;
}


/**
 * The next 3 utility functions are lifted from IOR.
 */
/******************************************************************************/
/*
 * Extract key/value pair from hint string.
 */

void
ExtractHint(char * settingVal,
            char * valueVal,
            char * hintString)
{
    char * settingPtr,
         * valuePtr,
         * tmpPtr1,
         * tmpPtr2;

    settingPtr = (char *)strtok(hintString, "=");
    valuePtr = (char *)strtok(NULL, " \t\r\n");
    tmpPtr1 = settingPtr;
    tmpPtr2 = (char *)strstr(settingPtr, "MPIIO_HINT__");
    if (tmpPtr1 == tmpPtr2) {
        settingPtr += strlen("MPIIO_HINT__");
    }
    strcpy(settingVal, settingPtr);
    strcpy(valueVal, valuePtr);
} /* ExtractHint() */


/******************************************************************************/
/*
 * Set hints for MPIIO, HDF5, or NCMPI.
 */
#define MAX_HINT_STR 1024
void
SetHints(MPI_Info * mpiHints, char * hintsFileName)
{
    char           hintString[MAX_HINT_STR],
                   settingVal[MAX_HINT_STR],
                   valueVal[MAX_HINT_STR];
    extern char ** environ;
    int            i;
    FILE         * fd;

    /*
     * This routine checks for hints from the environment and/or from the
     * hints files.  The hints are of the form:
     * 'MPIIO_HINT__<hint>=<value>', <hint> is the full name of the hint
     * to be set, and <value> is the hint value.
     * E.g., 'setenv MPIIO_HINT__panfs_concurrent_write 1'
     * or 'MPIIO_HINT__panfs_concurrent_write=1' in the hints file.
     */
    MPI_Info_create(mpiHints);

    /* get hints from environment */
    for (i = 0; environ[i] != NULL; i++) {
        /* if this is an MPIIO_HINT, pass the hint to the info object */
        if (strncmp(environ[i], "MPIIO_HINT", strlen("MPIIO_HINT")) == 0) {
            strcpy(hintString, environ[i]);
            ExtractHint(settingVal, valueVal, hintString);
            MPI_Info_set(*mpiHints, settingVal, valueVal);
        }
    }

    /* get hints from hints file */
    if (strcmp(hintsFileName, "") != 0) {

        /* open the hint file */
        fd = fopen(hintsFileName, "r");
        if (fd == NULL) {
            log_error(netcdf_debug_level, "cannot open hints file");
        } else {
            /* iterate over hints file */
            while(fgets(hintString, MAX_HINT_STR, fd) != NULL) {
                if (strncmp(hintString, "MPIIO_HINT", strlen("MPIIO_HINT")) == 0) {
                    ExtractHint(settingVal, valueVal, hintString);
                    MPI_Info_set(*mpiHints, settingVal, valueVal);
                }
            }
            /* close the hints files */
            if (fclose(fd) != 0) log_error(netcdf_debug_level, "cannot close hints file");
        }
    }
} /* SetHints() */


/******************************************************************************/
/*
 * Show all hints (key/value pairs) in an MPI_Info object.
 */

void ShowHints(MPI_Info * mpiHints)
{
    char key[MPI_MAX_INFO_VAL],
         value[MPI_MAX_INFO_VAL];
    int  flag,
         i,
         nkeys;

    MPI_Info_get_nkeys(*mpiHints, &nkeys);

    for (i = 0; i < nkeys; i++) {
        MPI_Info_get_nthkey(*mpiHints, i, key);
        MPI_Info_get(*mpiHints, key, MPI_MAX_INFO_VAL-1, value, &flag);
        log_debug(netcdf_debug_level, "mpiHint[%d]: %s = %s", i, key, value);
    }
} /* ShowHints() */

/* -------------------- SERVER-SIDE STUBS ---------- */

/**
 * @brief Create a netcdf dataset.
 *
 * This is the server-side stub that uses nc_create to create
 * a netcdf dataset.  It returns the
 */
int nc_create_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nc_create_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;
    int rank, np;
    const char *path = args->path;
    const int cmode = args->cmode;
    size_t chunksizehint = args->chunksizehint;
    int ncid=-1;
    nc_create_res res;  /* this is what we send back to the client */
    MPI_Info mpiHints = MPI_INFO_NULL;

    log_level debug_level = netcdf_debug_level;

    log_debug(debug_level, "calling nc_create(%s, %d, %d)", path, cmode, ncid);

    create_time=MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    SetHints(&mpiHints, "");
    ShowHints(&mpiHints);

    /* call the netcdf function */
    //rc = nc__create(path, cmode, initialsz, &chunksizehint, &ncid);

    double callTime;
    Start_Timer(callTime);
    rc = ncmpi_create(MPI_COMM_WORLD, path, cmode, mpiHints, &ncid);
    if (rc != NC_NOERR) {
        log_error(debug_level, "Could not create file(%s): %s", path, ncmpi_strerror(rc));
        goto cleanup;
    }
    Stop_Timer("create", callTime);

    /* set the result structure */
    memset(&res, 0, sizeof(res));
    res.chunksizehint = chunksizehint;
    res.ncid = ncid;

    if (args->num_participants > 0) {
        set_participant_count(ncid, caller, args->write_type, args->num_participants);
    } else {
        add_participant_for_file(ncid, caller, args->write_type);
    }

cleanup:
    log_debug(netcdf_debug_level, "exit (ncid=%d)", ncid);

    /* send the ncipd and return code back to client */
    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    return rc;
}

/**
 * Get the metadata associated a complete group hierarchy.
 *
 * A group consists of an id, name, dimensions (for all group vars),
 * and attributes.  This function gathers each of these attributes
 * through the nc_inq_* commands, then recursively gets the attributes
 * for each of the subgroups.
 */
int get_group_metadata(
        const int ncid,
        struct nc_group *group)
{
    int rc = NC_NOERR, i, j, natts;
    log_level debug_level = netcdf_debug_level;

    int *dimids=NULL;
    int *varids=NULL;
    int *groupids = NULL;

    /* set the ncid */
    group->ncid = ncid;


#ifdef USING_NETCDF4

    /* get the namelen for the name of this group */
    if ((rc = nc_inq_grpname_len(ncid, &namelen)) != NC_NOERR) {
        log_error(debug_level, "Could not get namelen: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

    /* get the name */
    group->name = (char *)calloc(namelen+1, sizeof(char));
    if ((rc = nc_inq_grpname_full(ncid, &namelen, group->name)) != NC_NOERR) {
        log_error(debug_level, "Could not get groupname: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

    /* get parent ncid */
    rc = nc_inq_grp_parent(ncid, &group->parent_ncid);
    if (rc == NC_ENOTNC4) {
        group->parent_ncid = -1;
        rc = NC_NOERR;
    }
    else if (rc != NC_NOERR) {
        log_error(debug_level, "Could not get parent ncid: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

#else

    /* Without group support (in NC4), there is only one group */
    group->name = strdup("/");
    group->parent_ncid = -1;


    /* get the unlimdimid */
    if ((rc = ncmpi_inq_unlimdim(ncid, &group->unlimdimid)) != NC_NOERR) {
        log_error(debug_level, "Error getting unlimdimid: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

#endif

    /* ----- Get Dimension Metadata for the group ----- */
    dimids = (int *)calloc(NC_MAX_DIMS, sizeof(int));

#ifdef USING_NETCDF4
    /* get the dimids for this group */
    if ((rc = nc_inq_dimids(ncid, (int *)&group->dims.dims_len, dimids, 0)) != NC_NOERR) {
        log_error(debug_level, "Error getting dimids: %s", ncmpi_strerror(rc));
        goto cleanup;
    }
#else
    /* get the number of dimids for this dataset */
    if ((rc = ncmpi_inq_ndims(ncid, (int *)&group->dims.dims_len)) != NC_NOERR) {
        log_error(debug_level, "Error getting dimids: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

    /* in netcdf3, dimids range from 0 -- (ndims-1) */
    for (i=0; i<group->dims.dims_len; i++) {
        dimids[i] = i;
    }
#endif

    log_debug(debug_level, "got ndims=%d (ncid=%d)", group->dims.dims_len, ncid);

    /* get metadata for each dimension */
    group->dims.dims_val = (nc_dim *)calloc(group->dims.dims_len, sizeof(nc_dim));
    for (i=0; i<group->dims.dims_len; i++) {

        /* set current dimension */
        struct nc_dim *dim = &group->dims.dims_val[i];
        dim->dimid = dimids[i];

        /* allocate space for name (assume largest) */
        dim->name = (char *)calloc(NC_MAX_NAME, sizeof(char));

        /* get dimension metadata */
        rc = ncmpi_inq_dim(ncid, dim->dimid, dim->name, (MPI_Offset*)&dim->len);
        if (rc != NC_NOERR) {
            log_error(debug_level, "Could not get dimid=%d metadata: %s", dimids[i], ncmpi_strerror(rc));
            goto cleanup;
        }
    }

    /* free the dimids, we will use them again for variables */
    free(dimids);
    dimids = NULL;

    /* ------- Global attribute metadata -------- */
    ncmpi_inq_natts(ncid, &natts);
    group->atts.atts_len = natts;
    group->atts.atts_val = (nc_att *)calloc(natts, sizeof(nc_att));

    for (i=0; i< natts; i++) {
        nc_att *att = &group->atts.atts_val[i];
        att->name = (char *)calloc(NC_MAX_NAME, sizeof(char));
        rc = ncmpi_inq_attname(ncid, NC_GLOBAL, i, att->name);
        if (rc != NC_NOERR) {
            log_error(debug_level, "Could not get global attributes for ncid=%d", ncid);
            goto cleanup;
        }
        rc = ncmpi_inq_att(ncid, NC_GLOBAL, att->name, (nc_type *)&att->xtype, (MPI_Offset*)&att->len);
        if (rc != NC_NOERR) {
            log_error(debug_level, "Could not get global attributes for ncid=%d", ncid);
            goto cleanup;
        }
        log_debug(netcdf_debug_level, "attr: name(%s) type(%d) len(%d)", att->name, att->xtype, att->len);
    }

    log_debug(debug_level, "Got attribute metadata");

    /* ----- Get Variable Metadata ----- */
    varids = (int *)calloc(NC_MAX_VARS, sizeof(int));

#ifdef USING_NETCDF4
    /* get varids for this group */
    if ((rc = nc_inq_varids(ncid, (int *)&group->vars.vars_len, varids)) != NC_NOERR) {
        log_error(debug_level, "Could not get varids: %s", ncmpi_strerror(rc));
        goto cleanup;
    }
#else
    /* get varids for this group */
    if ((rc = ncmpi_inq_nvars(ncid, (int *)&group->vars.vars_len)) != NC_NOERR) {
        log_error(debug_level, "Could not get varids: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

    /* in netcdf3, varids ranged from 0 to (nids-1) */
    for (i=0; i<group->vars.vars_len; i++) {
        varids[i] = i;
    }
#endif

    /* allocate space for variables */
    group->vars.vars_val = (nc_var *)calloc(group->vars.vars_len, sizeof(nc_var));

    log_debug(debug_level, "Got varids (%d vars)", (int)group->vars.vars_len);

    /* get metadata about each variable */
    for (i=0; i<group->vars.vars_len; i++) {

        /* Variable structure for this ID */
        struct nc_var *var = &group->vars.vars_val[i];
        var->varid = varids[i];

        /* First we need to allocate space for name and the dimids */
        var->name = (char *)calloc(NC_MAX_NAME, sizeof(char));
        var->dimids.dimids_val = (int *)calloc(NC_MAX_VAR_DIMS, sizeof(int));

        /* Get metadata for the variable */
        rc = ncmpi_inq_var(ncid, var->varid, var->name, (nc_type *)&var->xtype,
                (int *)&var->dimids.dimids_len, var->dimids.dimids_val,
                (int *)&var->atts.atts_len);
        if (rc != NC_NOERR) {
            log_error(debug_level, "Could not call nc_inq_var(varid=%d): %s",
                    var->varid, ncmpi_strerror(rc));
            goto cleanup;
        }

        log_debug(debug_level, "Got info for var[%d], var->varid=%d, var->name=%s", i, var->varid, var->name);

        /* Get the attributes */
        var->atts.atts_val = (nc_att *)calloc(var->atts.atts_len, sizeof(nc_att));
        for (j = 0; j<var->atts.atts_len; j++) {

            struct nc_att *att = &var->atts.atts_val[j];

            /* get attribute name */
            att->name = (char *)calloc(NC_MAX_NAME, sizeof(char));
            if ((rc = ncmpi_inq_attname(ncid, var->varid, j, att->name)) != NC_NOERR) {
                log_error(debug_level, "Could not get attribute name: %s", ncmpi_strerror(rc));
                goto cleanup;
            }

            /* get attribute info */
            if ((rc = ncmpi_inq_att(ncid, var->varid, att->name, (nc_type *)&att->xtype, (MPI_Offset*)&att->len)) != NC_NOERR) {
                log_error(debug_level, "Could not get attribute metadata: %s", ncmpi_strerror(rc));
                goto cleanup;
            }
        }
    }

    log_debug(debug_level, "Got variable metadata");

#ifdef USING_NETCDF4

    /* Get number of subgroups -- WHY NOT HAVE NC_INQ_NGRPS? */
    rc = nc_inq_grps(ncid, (int *)&group->groups.groups_len, NULL);
    if (rc == NC_ENOTNC4) {
        group->groups.groups_len = 0;
        rc = NC_NOERR;
    }
    else if (rc != NC_NOERR) {
        log_error(debug_level, "Unable to get grouplen: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

    log_debug(debug_level, "Got number of subgroups (%d)", (int)group->groups.groups_len);

    if (group->groups.groups_len > 0) {
        /* Get group IDs */
        groupids = (int *)calloc(group->groups.groups_len, sizeof(int));
        if ((rc = nc_inq_grps(ncid, NULL, groupids)) != NC_NOERR) {
            log_error(debug_level, "Unable to get grouplen: %s", ncmpi_strerror(rc));
            goto cleanup;
        }

        /* Recursively get subgroup metadata */
        group->groups.groups_val = (nc_group *)calloc(group->groups.groups_len, sizeof(nc_group));
        for (i=0; i<group->groups.groups_len; i++) {
            if ((rc = get_group_metadata(groupids[i], &group->groups.groups_val[i])) != NC_NOERR) {
                log_error(debug_level, "Unable to get subgroup metadata: %s", ncmpi_strerror(rc));
                goto cleanup;
            }
        }

        log_debug(debug_level, "Got all subgroup metadata");
    }
#endif


cleanup:
    if (dimids != NULL) free(dimids);
    if (varids != NULL) free(varids);
    if (groupids != NULL) free(groupids);

    return rc;
}


/**
 * @brief Open a netcdf dataset.
 *
 * Open a netcdf dataset and return all associated metadata.
 */
int nc_open_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nc_open_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;
    const char *path = args->path;
    const int omode = args->mode;
    size_t chunksizehint = args->chunksizehint;
    int ncid=-1;
    nc_open_res res;  /* this is what we send back to the client */
    MPI_Info mpiHints = MPI_INFO_NULL;

    memset(&res, 0, sizeof(res));

    log_level debug_level = netcdf_debug_level;

    log_debug(debug_level, "calling nc__open(%s, %d)", path, omode);

    /* call the netcdf function */
    /*
    if ((rc = nc__open(path, omode, &chunksizehint, &ncid)) != NC_NOERR) {
        log_error(debug_level, "Error opening file \"%s\": %s", path, ncmpi_strerror(rc));
        goto cleanup;
    }
    */

    open_time=MPI_Wtime();

    SetHints(&mpiHints, "");
    ShowHints(&mpiHints);

    double callTime;
    Start_Timer(callTime);
    if ((rc = ncmpi_open(MPI_COMM_WORLD, path, omode, mpiHints, &ncid)) != NC_NOERR) {
        log_error(debug_level, "Error opening file \"%s\": %s", path, ncmpi_strerror(rc));
        res.root_group.ncid = -1;
        goto cleanup;
    }
    Stop_Timer("open", callTime);

    log_debug(debug_level, "dataset open (ncid=%d)", ncid);

    /* set the result structure */
    res.chunksizehint = chunksizehint;
    res.root_group.ncid = ncid;

    //root_group = (struct nc_group *)calloc(1, sizeof(struct nc_group));

    rc = get_group_metadata(ncid, &res.root_group);
    if (rc != NC_NOERR) {
        log_error(debug_level, "Could not get root_group metadata");
        goto cleanup;
    }

    if (args->num_participants > 0) {
        set_participant_count(ncid, caller, args->write_type, args->num_participants);
    } else {
        add_participant_for_file(ncid, caller, args->write_type);
    }

cleanup:
    log_debug(netcdf_debug_level, "exit (ncid=%d)", ncid);

    /* send the ncid and return code back to client */
    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    return rc;
}

/**
 * @brief Define a dimension of a dataset.
 *
 * This is the server-side stub that uses nc_create to create
 * a netcdf dataset.  It returns the
 */
int nc_def_dim_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nc_def_dim_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;
    const int ncid = args->ncid;
    const char *name = args->name;
    const size_t len = args->len;
    int dimid; /* result */

    log_debug(netcdf_debug_level, "enter");

    double callTime;
    Start_Timer(callTime);
    /* call actual netcdf function */
    rc = ncmpi_def_dim(ncid, name, len, &dimid);
    Stop_Timer("def dim", callTime);

    log_debug(netcdf_debug_level, "nc_def_dim: rc=%s", ncmpi_strerror(rc));

    /* send the ncipd and return code back to client */
    rc = nssi_send_result(caller, request_id, rc, &dimid, res_addr);

    log_debug(netcdf_debug_level, "finished nc_def_dim(ncid=%d, name=%s, len=%d, dimid=%d)",
            ncid, name, len, dimid);

    log_debug(netcdf_debug_level, "exit (dimid=%d)", dimid);

    return rc;
}

/**
 * @brief Define a variable of a dataset.
 *
 * This is the server-side stub that uses nc_create to create
 * a netcdf dataset.  It returns the
 */
int nc_def_var_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nc_def_var_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    log_debug(netcdf_debug_level, "enter");

    log_debug(netcdf_debug_level, "calling nc_def_var");

    int rc = NC_NOERR;
    const int ncid = args->ncid;
    const char *name = args->name;
    const nc_type xtype = (nc_type)args->xtype;
    const int ndims = args->dimids.dimids_len;
    const int *dimids = args->dimids.dimids_val;
    int varid;  /* result */

    log_debug(netcdf_debug_level, "calling nc_def_var(ncid=%d, name=\"%s\", xtype=%d, "
            "ndims=%d, dimids[0]=%d)", ncid, name, xtype, ndims, dimids[0]);


    double callTime;
    Start_Timer(callTime);
    /* call real netcdf function */
    rc = ncmpi_def_var(ncid, name, xtype, ndims, dimids, &varid);
    if (rc != NC_NOERR) {
        log_error(netcdf_debug_level, "%s", ncmpi_strerror(rc));
    }
    Stop_Timer("def var", callTime);

    /* send the ncipd and return code back to client */
    rc = nssi_send_result(caller, request_id, rc, &varid, res_addr);

    log_debug(netcdf_debug_level, "exit (varid=%d)", varid);

    return rc;
}

int ncmpi_inq_type(
        const int ncid,
        const nc_type xtype,
        char *name,
        MPI_Offset *sizep)
{
    int rc = NC_NOERR;

    log_debug(netcdf_debug_level, "enter");

    switch (xtype) {

    case NC_BYTE:
        *sizep = 1;
        if (name) strcpy(name, "NC_BYTE");
        break;

    case NC_CHAR:
        *sizep = sizeof(char);
        if (name) strcpy(name, "NC_CHAR");
        break;

    case NC_SHORT:
        *sizep = sizeof(short);
        if (name) strcpy(name, "NC_SHORT");
        break;

    case NC_INT:
        *sizep = sizeof(int);
        if (name) strcpy(name, "NC_INT");
        break;

    case NC_FLOAT:
        *sizep = sizeof(float);
        if (name) strcpy(name, "NC_FLOAT");
        break;

    case NC_DOUBLE:
        *sizep = sizeof(double);
        if (name) strcpy(name, "NC_DOUBLE");
        break;

    default:
        rc = NC_EBADTYPE;
        break;
    }

    log_debug(netcdf_debug_level, "exit");

    return rc;
}

int nc_get_att_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nc_get_att_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;
    log_level debug_level = netcdf_debug_level;

    const int ncid = args->ncid;
    const int varid = args->varid;
    const char *name = args->name;

    nc_type xtype;
    MPI_Offset len, nbytes, typesize;
    void *attbuf = NULL;

    log_debug(netcdf_debug_level, "enter");

    /* get the len and type of the attribute buffer */
    rc = ncmpi_inq_att(ncid, varid, name, &xtype, &len);
    if (rc != NC_NOERR) {
        log_error(debug_level, "Unable to inq_att: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

    /* get the size of the datatype */
    rc = ncmpi_inq_type(ncid, xtype, NULL, &typesize);
    if (rc != NC_NOERR) {
        log_error(debug_level, "Unable to inq_type: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

    nbytes = len * typesize;
    attbuf = calloc(nbytes, sizeof(char));

    switch (xtype) {

    case NC_CHAR:
        rc = ncmpi_get_att_text(ncid, varid, name, (char *)attbuf);
        break;

    case NC_SHORT:
        rc = ncmpi_get_att_short(ncid, varid, name, (short *)attbuf);
        break;

    case NC_INT:
        rc = ncmpi_get_att_int(ncid, varid, name, (int *)attbuf);
        break;

    case NC_FLOAT:
        rc = ncmpi_get_att_float(ncid, varid, name, (float *)attbuf);
        break;

    case NC_DOUBLE:
        rc = ncmpi_get_att_double(ncid, varid, name, (double *)attbuf);
        break;

    default:
        rc = NC_EIOMISMATCH;
        break;

    }

    if (rc != NC_NOERR) {
        log_error(debug_level, "%s", ncmpi_strerror(rc));
        goto cleanup;
    }


    /* send the buf to the client */
    rc = nssi_put_data(caller, attbuf, nbytes, data_addr, -1);
    if (rc != NSSI_OK) {
        log_error(debug_level, "%s", nssi_err_str(rc));
        goto cleanup;
    }

cleanup:

    /* send result back to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    if (attbuf) free(attbuf);

    log_debug(netcdf_debug_level, "exit");

    return rc;
}

int nc_put_att_stub(
        const unsigned long request_id,
                const NNTI_peer_t *caller,
                const nc_put_att_args *args,
                const NNTI_buffer_t *data_addr,
                const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;

    const int ncid = args->ncid;
    const int varid = args->varid;
    const char *name = args->name;
    const nc_type xtype = (nc_type)args->xtype;  /* one of the 6 netcdf xtypes */
    const arg_type atype = args->atype;
    const ssize_t len = args->data.data_len;
    const char *data = args->data.data_val;

    log_debug(netcdf_debug_level, "enter");


    /* atype (arg type) identifies the function to call */
    switch(atype) {

        case NC_ARG_TEXT:
        {
            log_debug(netcdf_debug_level, "calling nc_put_att_text(ncid=%d, varid=%d, "
                    "name=%s, len=%d, data=%s)", ncid, varid, name, len, data);

            rc = ncmpi_put_att_text(ncid, varid, name,
                    len/sizeof(char), data);
            break;
        }

        case NC_ARG_UCHAR:
        {
            rc = ncmpi_put_att_uchar(ncid, varid, name, xtype,
                    len/sizeof(unsigned char), (const unsigned char *)data);
            break;
        }

        case NC_ARG_SCHAR:
        {
            rc = ncmpi_put_att_schar(ncid, varid, name, xtype,
                    len/sizeof(signed char), (const signed char *)data);
            break;
        }

        case NC_ARG_SHORT:
        {
            rc = ncmpi_put_att_short(ncid, varid, name, xtype,
                    len/sizeof(short), (const short *)data);
            break;
        }

        case NC_ARG_INT:
        {
            rc = ncmpi_put_att_int(ncid, varid, name, xtype,
                    len/sizeof(int), (const int *)data);
            break;
        }

        case NC_ARG_LONG:
        {
            rc = ncmpi_put_att_long(ncid, varid, name, xtype,
                    len/sizeof(long), (const long *)data);
            break;
        }

        case NC_ARG_FLOAT:
        {
            rc = ncmpi_put_att_float(ncid, varid, name, xtype,
                    len/sizeof(float), (const float *)data);
            break;
        }

        case NC_ARG_DOUBLE:
        {
            rc = ncmpi_put_att_double(ncid, varid, name, xtype,
                    len/sizeof(double), (const double *)data);
            break;
        }

        case NC_ARG_UBYTE:
        {
#if HAVE_NC_PUT_ATT_UBYTE
            rc = ncmpi_put_att_ubyte(ncid, varid, name, xtype,
                    len/sizeof(unsigned char), (const unsigned char *)data);
#else
            rc = NC_ENOTSUPP;  /* error added to handle lack of support */
#endif
            break;
        }

        case NC_ARG_USHORT:
        {
#if HAVE_NC_PUT_ATT_USHORT
            rc = ncmpi_put_att_ushort(ncid, varid, name, xtype,
                    len/sizeof(unsigned short), (const unsigned short *)data);
#else
            rc = NC_ENOTSUPP;
#endif
            break;
        }

        case NC_ARG_UINT:
        {
#if HAVE_NC_PUT_ATT_UINT
            rc = ncmpi_put_att_uint(ncid, varid, name, xtype,
                    len/sizeof(unsigned int), (const unsigned int *)data);
#else
            rc = NC_ENOTSUPP;
#endif
            break;
        }

        case NC_ARG_LONGLONG:
        {
#if HAVE_NC_PUT_ATT_LONGLONG
            rc = ncmpi_put_att_longlong(ncid, varid, name, xtype,
                    len/sizeof(long long), (const long long *)data);
#else
            rc = NC_ENOTSUPP;
#endif
            break;
        }

        case NC_ARG_ULONGLONG:
        {
#if HAVE_NC_PUT_ATT_ULONGLONG
            rc = ncmpi_put_att_ulonglong(ncid, varid, name, xtype,
                    len/sizeof(unsigned long long), (const unsigned long long *)data);
#else
            rc = NC_ENOTSUPP;
#endif
            break;
        }

        case NC_ARG_VOID:
        {
            //rc = ncmpi_put_att(ncid, varid, name, xtype, len, data);
            //break;
        }

        default:
            rc = NC_EBADTYPE;
    }


    /* send result back to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    log_debug(netcdf_debug_level, "exit");

    return rc;
}


/**
 * Write a variable to the underlying netcdf file.
 *
 * This implementation will pull raw data from the client
 * using the data.
 *
 * TODO: Put data in a cache before calling the actual netcdf function.
 *
 */
int nc_put_vars_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nc_put_vars_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;
    int i, ndims;

    const int ncid = args->ncid;
    const int varid = args->varid;
    const int atype = args->atype;
    const nc_type buftype = (nc_type)args->buftype;
    const size_t count    = args->element_count;
    const size_t len      = args->len;
    const nc_size_t *startp = args->start.start_val;
    const nc_size_t *countp = args->count.count_val;
    const nc_size_t *stridep = args->stride.stride_val;
    MPI_Offset *count_copy = NULL;
    MPI_Offset *start_copy = NULL;
    MPI_Offset *stride_copy = NULL;
    int *dimids = NULL;

    log_level debug_level = netcdf_debug_level;

    MPI_Datatype datatype;

    int ccount = 1;

    void *buf=NULL;


    log_debug(netcdf_debug_level, "enter");

    log_debug(debug_level, "calling nc_put_vars_stub(ncid=%d, varid=%d, atype=%d, len=%d)", ncid, varid, atype, len);

    double callTime;

    Start_Timer(callTime);
    buf = malloc(len);
    Stop_Timer("buf malloc", callTime);

    Start_Timer(callTime);

    rc = ncmpi_inq_varndims(ncid, varid, &ndims);
    if (rc != NC_NOERR) {
        log_error(debug_level, "Could not get ndims: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

    if (ndims) {
        size_t nbytes = 0;

        /* Treat each arg type differently */
        switch (buftype) {
            case NC_BYTE:
                datatype = MPI_BYTE;
                log_debug(debug_level, "using MPI_BYTE (count=%d ; len=%d)", count, len);
                break;
            case NC_CHAR:
                datatype = MPI_CHAR;
                log_debug(debug_level, "using MPI_CHAR (count=%d ; len=%d)", count, len);
                break;

            case NC_SHORT:
                datatype = MPI_SHORT;
                log_debug(debug_level, "using MPI_SHORT (count=%d ; len=%d)", count, len);
                break;

            case NC_INT:
                datatype = MPI_INT;
                log_debug(debug_level, "using MPI_INT (count=%d ; len=%d)", count, len);
                break;

            case NC_FLOAT:
                datatype = MPI_FLOAT;
                log_debug(debug_level, "using MPI_FLOAT (count=%d ; len=%d)", count, len);
                break;

            case NC_DOUBLE:
                datatype = MPI_DOUBLE;
                log_debug(debug_level, "using MPI_DOUBLE (count=%d ; len=%d)", count, len);
                break;
            default:
                log_error(debug_level, "Operation=%d not supported", buftype);
                rc = NC_ENOTSUPP;
                goto cleanup;
        }
        start_copy = (MPI_Offset *)calloc(ndims, sizeof(MPI_Offset));
        count_copy = (MPI_Offset *)calloc(ndims, sizeof(MPI_Offset));
        stride_copy = (MPI_Offset *)calloc(ndims, sizeof(MPI_Offset));

        /* copy values */
        if (startp) {
            copy(startp, startp+ndims, start_copy);
        }
        else {
            /* default to 0,0,...,0 */
            fill(start_copy, start_copy+ndims, 0);
        }

        /* set counts */
        if (countp) {
            copy(countp, countp+ndims, count_copy);
        }
        else {
            /* defaults to dimlen[0], dimlen[1], ... dimlen[ndims-1] */

            dimids = (int *) calloc(ndims, sizeof(int));

            /* get the diminsion IDs for this variable */
            rc = ncmpi_inq_vardimid(ncid, varid, dimids);
            if (rc != NC_NOERR) {
                log_error(debug_level, "could not get dimids");
                goto cleanup;
            }

            /* set count[i] to dimlen[i] */
            for (i = 0; i<ndims; i++) {
                MPI_Offset dimlen;

                rc = ncmpi_inq_dimlen(ncid, dimids[i], &dimlen);
                if (rc != NC_NOERR) {
                    log_error(debug_level, "could not get dimlen");
                    goto cleanup;
                }

                log_debug(debug_level, "dimid[%d] = %d", dimids[i], dimlen);

                count_copy[i] = dimlen;
            }
        }

        if (stridep) {
            copy(stridep, stridep+ndims, stride_copy);
        }
        else {
            /* defaults to 1,1,1... */
            fill(stride_copy, stride_copy+ndims, 1);
        }

        fprintf(logger_get_file(), "start_copy = (");
        for (i=0;i<ndims; i++) {
            fprintf(logger_get_file(), "%s%d", !i ? "" : ", ", start_copy[i]);
        }
        fprintf(logger_get_file(), ")\n");
        fprintf(logger_get_file(), "count_copy = (");
        for (i=0;i<ndims; i++) {
            fprintf(logger_get_file(), "%s%d", !i ? "" : ", ", count_copy[i]);
        }
        fprintf(logger_get_file(), ")\n");
        fprintf(logger_get_file(), "stride_copy = (");
        for (i=0;i<ndims; i++) {
            fprintf(logger_get_file(), "%s%d", !i ? "" : ", ", stride_copy[i]);
        }
        fprintf(logger_get_file(), ")\n");

        ccount = 1;
        for (i=0;i<ndims; i++) {
            ccount *= count_copy[i];  /* product of the dimlens */
        }

        /* send error if our calculated count does not equal the count sent by the client */
        if (ccount != count) {
            rc = NC_EIOMISMATCH;
            log_error(debug_level, "%s", ncmpi_strerror(rc));
            goto cleanup;
        }

        /* TODO: put the rest of this op into an async job */

        log_debug(debug_level, "Fetch %d vals (%lu bytes) from client", count, len);
        /* Fetch the data from the client */
        rc = nssi_get_data(caller, buf, len, data_addr);
        if (rc != NSSI_OK) {
            log_error(debug_level, "Could not fetch var data from client");
            goto cleanup;
        }

        if (use_aggregation(ncid)) {
            log_debug(netcdf_debug_level, "using agg");
            aggregation_chunk_details_t *chunk;
            chunk = new aggregation_chunk_details_t;
            chunk->ncid = ncid;
            chunk->varid = varid;
            chunk->ndims = ndims;
            chunk->buf = buf;
            chunk->atype = atype;
            chunk->len   = len;
            chunk->datatype      = datatype;
            chunk->num_elements  = count;
            chunk->start  = start_copy;
            chunk->count  = count_copy;
            chunk->stride = stride_copy;
            add_participant_chunk(caller, chunk);

            rc = NSSI_OK;

            if (aggregate_data_ready_to_write(ncid, varid)) {
                log_debug(netcdf_debug_level, "ready to agg");
                int chunk_count=0;

                double aggTime;
                Start_Timer(aggTime);
                while(try_aggregation(ncid, varid) == TRUE);
                Stop_Timer("agg", aggTime);

                aggregation_chunk_details_t **chunks = get_chunks(ncid, varid, &chunk_count);
                if (chunk_count == 0) {
                    log_error(netcdf_debug_level, "chunk_count is 0.  how can that be?");
                    goto cleanup;
                }

                if (use_collective(ncid)) {
                    log_debug(netcdf_debug_level, "using collective");
                    double collPutTime;
                    Start_Timer(collPutTime);
                    rc = do_put_vars_all(ncid, varid, chunks, chunk_count);
                    Stop_Timer("coll put", collPutTime);
                    if (rc != NC_NOERR) {
                        log_error(debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
                        free(chunks);
                        goto cleanup;
                    }
                } else {
                    log_debug(netcdf_debug_level, "using independent");
                    double indepPutTime;
                    Start_Timer(indepPutTime);
                    rc = do_put_vars(chunks, chunk_count);
                    Stop_Timer("indep put", indepPutTime);
                    if (rc != NC_NOERR) {
                        log_error(debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
                        free(chunks);
                        goto cleanup;
                    }
                }

                free(chunks);

                cleanup_aggregation_chunks(ncid, varid);
            } else {
                log_debug(netcdf_debug_level, "not ready to agg");
            }
        } else if (use_caching(ncid)) {
            log_debug(netcdf_debug_level, "using caching");
            aggregation_chunk_details_t *chunk;
            chunk = new aggregation_chunk_details_t;
            chunk->ncid = ncid;
            chunk->varid = varid;
            chunk->ndims = ndims;
            chunk->buf = buf;
            chunk->atype = atype;
            chunk->len   = len;
            chunk->datatype      = datatype;
            chunk->num_elements  = count;
            chunk->start  = start_copy;
            chunk->count  = count_copy;
            chunk->stride = stride_copy;
            add_participant_chunk(caller, chunk);

            rc = NSSI_OK;

            if (cache_data_ready_to_write(ncid, varid)) {
                log_debug(netcdf_debug_level, "ready to flush cache");
                int chunk_count=0;

                aggregation_chunk_details_t **chunks = get_chunks(ncid, varid, &chunk_count);
                if (chunk_count == 0) {
                    log_error(netcdf_debug_level, "chunk_count is 0.  how can that be?");
                    goto cleanup;
                }

                if (use_collective(ncid)) {
                    log_debug(netcdf_debug_level, "using collective");
                    double collPutTime;
                    Start_Timer(collPutTime);
                    rc = do_put_vars_all(ncid, varid, chunks, chunk_count);
                    Stop_Timer("coll put", collPutTime);
                    if (rc != NC_NOERR) {
                        log_error(debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
                        free(chunks);
                        goto cleanup;
                    }
                } else {
                    log_debug(netcdf_debug_level, "using independent");
                    double indepPutTime;
                    Start_Timer(indepPutTime);
                    rc = do_put_vars(chunks, chunk_count);
                    Stop_Timer("indep put", indepPutTime);
                    if (rc != NC_NOERR) {
                        log_error(debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
                        free(chunks);
                        goto cleanup;
                    }
                }

                free(chunks);

                cleanup_aggregation_chunks(ncid, varid);
            } else {
                log_debug(netcdf_debug_level, "not ready to flush cache");
            }
        } else {
            log_debug(netcdf_debug_level, "using direct");
            rc = do_put_vars(ncid, varid, start_copy, count_copy, stride_copy,
                             buf, count, datatype);
            if (rc != NC_NOERR) {
                log_error(debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
                goto cleanup;
            }
        }
    }

cleanup:
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);
    if (dimids) free(dimids);

    /* TODO: add the job to the job queue */
    if ((!use_aggregation(ncid)) &&
        (!use_caching(ncid))) {
        if (count_copy) free(count_copy);
        if (start_copy) free(start_copy);
        if (stride_copy) free(stride_copy);
        free(buf);
    }

    Stop_Timer("put var", callTime);

    log_debug(netcdf_debug_level, "exit");

    return rc;
}

/**
 * Read a variable to the underlying netcdf file.
 *
 * This implementation will push raw data to the client
 * using the data address.
 *
 * TODO: Put data in a cache before calling the actual netcdf function.
 *
 */
int nc_get_vars_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nc_get_vars_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;
    int i, ndims;

    const int ncid = args->ncid;
    const int varid = args->varid;
    const int atype = args->atype;
    const nc_size_t *stridep = args->stride.stride_val;
    const nc_size_t *startp = args->start.start_val;
    const nc_size_t *countp = args->count.count_val;
    const ssize_t len = args->len;

    MPI_Offset *count_copy = NULL;
    MPI_Offset *start_copy = NULL;
    MPI_Offset *stride_copy = NULL;
    int *dimids = NULL;

    nc_type vartype;
    log_level debug_level = netcdf_debug_level;

    MPI_Datatype datatype;
    int count = 1, ccount=0;
    void *buf = malloc(len);


    log_debug(netcdf_debug_level, "enter");

    log_debug(debug_level, "calling nc_get_vars_stub(ncid=%d, varid=%d, atype=%d, len=%d)", ncid, varid, atype, len);

    double callTime;
    Start_Timer(callTime);
    if (use_aggregation(ncid)) {
        flush_aggregated_chunks(ncid, varid);
    } else if (use_caching(ncid)) {
        flush_cached_chunks(ncid, varid);
    }
    Stop_Timer("flush on sync", callTime);


    rc = ncmpi_inq_varndims(ncid, varid, &ndims);
    if (rc != NC_NOERR) {
        log_error(debug_level, "Could not get ndims: %s", ncmpi_strerror(rc));
        goto cleanup;
    }

    if (ndims) {
        size_t count = 1;
        size_t nbytes = 0;

        /* find out the right type for this variable */
        rc = ncmpi_inq_vartype(ncid, varid, &vartype);
        if (rc != NC_NOERR) {
            log_error(debug_level, "Unable to get variable type");
            goto cleanup;
        }

        /* Treat each arg type differently */
        switch (vartype) {

        case NC_BYTE:
            datatype = MPI_BYTE;
            count = len;
            log_debug(debug_level, "using MPI_BYTE");
            break;

        case NC_CHAR:
            datatype = MPI_CHAR;
            count = len/sizeof(char);
            log_debug(debug_level, "using MPI_CHAR");
            break;

        case NC_SHORT:
            datatype = MPI_SHORT;
            count = len/sizeof(short);
            log_debug(debug_level, "using MPI_SHORT");
            break;

        case NC_INT:
            datatype = MPI_INT;
            count = len/sizeof(int);
            log_debug(debug_level, "using MPI_INT");
            break;

        case NC_FLOAT:
            datatype = MPI_FLOAT;
            count = len/sizeof(float);
            log_debug(debug_level, "using MPI_FLOAT");
            break;

        case NC_DOUBLE:
            datatype = MPI_DOUBLE;
            count = len/sizeof(double);
            log_debug(debug_level, "using MPI_DOUBLE");
            break;

        default:
            log_error(debug_level, "Operation=%d not supported", vartype);
            rc = NC_ENOTSUPP;
            goto cleanup;
        }

        log_debug(debug_level, "Sending %d vals (%lu bytes) to client", count, len);

        start_copy = (MPI_Offset *)calloc(ndims, sizeof(MPI_Offset));
        count_copy = (MPI_Offset *)calloc(ndims, sizeof(MPI_Offset));
        stride_copy = (MPI_Offset *)calloc(ndims, sizeof(MPI_Offset));

        /* copy values */
        if (startp) {
            copy(startp, startp+ndims, start_copy);
        }
        else {
            /* default to 0,0,...,0 */
            fill(start_copy, start_copy+ndims, 0);
        }

        /* set counts */
        if (countp) {
            copy(countp, countp+ndims, count_copy);
        }
        else {
            /* defaults to dimlen[0], dimlen[1], ... dimlen[ndims-1] */

            dimids = (int *) calloc(ndims, sizeof(int));

            /* get the diminsion IDs for this variable */
            rc = ncmpi_inq_vardimid(ncid, varid, dimids);
            if (rc != NC_NOERR) {
                log_error(debug_level, "could not get dimids");
                goto cleanup;
            }

            /* set count[i] to dimlen[i] */
            for (i = 0; i<ndims; i++) {
                MPI_Offset dimlen;

                rc = ncmpi_inq_dimlen(ncid, dimids[i], &dimlen);
                if (rc != NC_NOERR) {
                    log_error(debug_level, "could not get dimlen");
                    goto cleanup;
                }

                log_debug(debug_level, "dimid[%d] = %d", dimids[i], dimlen);

                count_copy[i] = dimlen;
            }
        }

        if (stridep) {
            copy(stridep, stridep+ndims, stride_copy);
        }
        else {
            /* defaults to 1,1,1... */
            fill(stride_copy, stride_copy+ndims, 1);
        }

        ccount = 1;
        for (i=0;i<ndims; i++) {
            ccount *= count_copy[i];  /* product of the dimlens */
        }

        /* send error if ccount does not equal count */
        if (ccount != count) {
            rc = NC_EIOMISMATCH;
            goto cleanup;
        }

        /* TODO: capture the rest of this operation as an async job request */

        rc = ncmpi_get_vars(ncid, varid, start_copy, count_copy, stride_copy,
                buf, count, datatype);
        if (rc != NC_NOERR) {
            log_error(debug_level, "Failed to put variable: %s", ncmpi_strerror(rc));
            goto cleanup;
        }

        /* push the data to the client */
        rc = nssi_put_data(caller, buf, len, data_addr, -1);
        if (rc != NSSI_OK) {
            log_error(debug_level, "Could not put var data on client");
            goto cleanup;
        }
    }

cleanup:
    /* we always send a result */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    if (buf) free(buf);

    if (dimids) free(dimids);

    /* TODO: add the request to the queue of job requests */

    /* this stuff gets freed when the job is complete */
    if (count_copy) free(count_copy);
    if (start_copy) free(start_copy);
    if (stride_copy) free(stride_copy);

    log_debug(netcdf_debug_level, "exit");

    return rc;
}



/**
 * @brief Restart dataset definition.
 */
int nc_redef_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const int *ncidp,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;

    log_debug(netcdf_debug_level, "enter");

    log_debug(netcdf_debug_level, "calling nc_redef(%d)", *ncidp);

    flush_aggregated_chunks(*ncidp);

    double callTime;
    Start_Timer(callTime);
    /* call real netcdf function */
    rc = ncmpi_redef(*ncidp);
    Stop_Timer("redef", callTime);

    /* send the ncipd and return code back to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    log_debug(netcdf_debug_level, "exit");

    return rc;
}

/**
 * @brief Complete dataset definition.
 */
int nc_enddef_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const int *ncidp,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;

    log_debug(netcdf_debug_level, "enter");

    log_debug(netcdf_debug_level, "calling nc_enddef(%d)", *ncidp);

    enddef_time=MPI_Wtime();

    double callTime;
    Start_Timer(callTime);
    /* call real netcdf function */
    rc = ncmpi_enddef(*ncidp);
    Stop_Timer("enddef", callTime);

    /* send the ncipd and return code back to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    Start_Timer(callTime);

    log_debug(netcdf_debug_level, "exit");

    return rc;
}

int nc_close_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const int *ncidp,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;

    log_debug(netcdf_debug_level, "enter");

    log_debug(netcdf_debug_level, "calling nc_close(%d)", *ncidp);

    close1_time=MPI_Wtime();

    double callTime;
    Start_Timer(callTime);
    if (use_aggregation(*ncidp)) {
        flush_aggregated_chunks(*ncidp);
    } else if (use_caching(*ncidp)) {
        flush_cached_chunks(*ncidp);
    }
    Stop_Timer("flush on close", callTime);

    Start_Timer(callTime);
    /* call real netcdf function */
    rc = ncmpi_close(*ncidp);
    if (rc != NC_NOERR) {
        log_error(netcdf_debug_level, "%s", ncmpi_strerror(rc));
    }
    Stop_Timer("close", callTime);

    remove_participant_for_file(*ncidp, caller);

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    close2_time=MPI_Wtime();


    log_debug(LOG_ALL, "recorded times create(%10.8f) open(%10.8f) enddef(%10.8f) close1(%10.8f) close2(%10.8f)",
            create_time, open_time, enddef_time, close1_time, close2_time);
    if (create_time == 0.0) {
        log_debug(LOG_ALL, "elapsed times enddef-open(%10.8f) close1-enddef(%10.8f) close2-enddef(%10.8f) close1-open(%10.8f) close2-open(%10.8f)",
                enddef_time-open_time, close1_time-enddef_time, close2_time-enddef_time, close1_time-open_time, close2_time-open_time);
    } else {
        log_debug(LOG_ALL, "elapsed times enddef-create(%10.8f) close1-enddef(%10.8f) close2-enddef(%10.8f) close1-create(%10.8f) close2-create(%10.8f)",
                enddef_time-create_time, close1_time-enddef_time, close2_time-enddef_time, close1_time-create_time, close2_time-create_time);
    }

    log_debug(netcdf_debug_level, "exit");

    return rc;
}

int nc_sync_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const int *ncidp,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;

    log_debug(netcdf_debug_level, "enter");

    log_debug(netcdf_debug_level, "calling nc_sync(%d)", *ncidp);

    double callTime;
    Start_Timer(callTime);
    if (use_aggregation(*ncidp)) {
        flush_aggregated_chunks(*ncidp);
    } else if (use_caching(*ncidp)) {
        flush_cached_chunks(*ncidp);
    }
    Stop_Timer("flush on sync", callTime);

    Start_Timer(callTime);
    /* call real netcdf function */
    rc = ncmpi_sync(*ncidp);
    if (rc != NC_NOERR) {
        log_error(netcdf_debug_level, "%s", ncmpi_strerror(rc));
    }
    Stop_Timer("sync", callTime);

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    log_debug(netcdf_debug_level, "exit");

    return rc;
}

int nc_begin_indep_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const int *ncidp,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;

    log_debug(netcdf_debug_level, "enter");

    log_debug(netcdf_debug_level, "calling nc_begin_indep_data(%d)", *ncidp);

    double callTime;
    Start_Timer(callTime);
    rc = ncmpi_begin_indep_data(*ncidp);
    Stop_Timer("begin indep", callTime);

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    log_debug(netcdf_debug_level, "exit");

    return rc;
}

int nc_end_indep_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const int *ncidp,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;

    log_debug(netcdf_debug_level, "enter");

    log_debug(netcdf_debug_level, "calling nc_end_indep_data(%d)", *ncidp);

    double callTime;
    Start_Timer(callTime);
    if (use_aggregation(*ncidp)) {
        flush_aggregated_chunks(*ncidp);
    } else if (use_caching(*ncidp)) {
        flush_cached_chunks(*ncidp);
    }
    Stop_Timer("flush on end indep", callTime);

    Start_Timer(callTime);
    rc = ncmpi_end_indep_data(*ncidp);
    Stop_Timer("end indep", callTime);

    /* send result to client */
    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    log_debug(netcdf_debug_level, "exit");

    return rc;
}

/**
 * @brief Create a netcdf dataset.
 *
 * This is the server-side stub that uses nc_create to create
 * a netcdf dataset.  It returns the
 */
int nc_set_fill_stub(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nc_set_fill_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NC_NOERR;
    int rank, np;
    const int ncid     = args->ncid;
    const int new_mode = args->new_fill_mode;
    nc_set_fill_res res;  /* this is what we send back to the client */

    log_level debug_level = netcdf_debug_level;

    log_debug(netcdf_debug_level, "enter");

    rc=ncmpi_set_fill(ncid, new_mode, &res.old_fill_mode);

    rc = nssi_send_result(caller, request_id, rc, &res, res_addr);

    log_debug(netcdf_debug_level, "exit");

    return rc;
}

/* -------- END SERVER-SIDE STUBS -------------- */

int netcdf_server_init()
{
    /* register server stubs */
    NSSI_REGISTER_SERVER_STUB(NETCDF_CREATE_OP, nc_create_stub, nc_create_args, nc_create_res);
    NSSI_REGISTER_SERVER_STUB(NETCDF_OPEN_OP, nc_open_stub, nc_open_args, nc_open_res);
    NSSI_REGISTER_SERVER_STUB(NETCDF_DEF_DIM_OP, nc_def_dim_stub, nc_def_dim_args, int);
    NSSI_REGISTER_SERVER_STUB(NETCDF_DEF_VAR_OP, nc_def_var_stub, nc_def_var_args, int);
    NSSI_REGISTER_SERVER_STUB(NETCDF_GET_ATT_OP, nc_get_att_stub, nc_get_att_args, void);
    NSSI_REGISTER_SERVER_STUB(NETCDF_PUT_ATT_OP, nc_put_att_stub, nc_put_att_args, void);
    NSSI_REGISTER_SERVER_STUB(NETCDF_REDEF_OP, nc_redef_stub, int, void);
    NSSI_REGISTER_SERVER_STUB(NETCDF_ENDDEF_OP, nc_enddef_stub, int, void);
    NSSI_REGISTER_SERVER_STUB(NETCDF_PUT_VARS_OP, nc_put_vars_stub, nc_put_vars_args, void);
    NSSI_REGISTER_SERVER_STUB(NETCDF_GET_VARS_OP, nc_get_vars_stub, nc_get_vars_args, void);
    NSSI_REGISTER_SERVER_STUB(NETCDF_SYNC_OP, nc_sync_stub, int, void);
    NSSI_REGISTER_SERVER_STUB(NETCDF_CLOSE_OP, nc_close_stub, int, void);
    NSSI_REGISTER_SERVER_STUB(NETCDF_BEGIN_INDEP_OP, nc_begin_indep_stub, int, void);
    NSSI_REGISTER_SERVER_STUB(NETCDF_END_INDEP_OP, nc_end_indep_stub, int, void);
    NSSI_REGISTER_SERVER_STUB(NETCDF_SET_FILL_OP, nc_set_fill_stub, nc_set_fill_args, nc_set_fill_res);

    return 0;
}



/* -- private methods -- */
static void print_opts(
        FILE *fp,
        const struct netcdf_args *args_info,
        const char *prefix)
{
    fprintf(fp, "%s ------------ NETCDF SERVER ---------\n", prefix);

    fprintf(fp, "%s ------------  ARGUMENTS -----------\n", prefix);
    fprintf(fp, "%s \tserver-url      = %s\n", prefix, args_info->server_url.c_str());
    fprintf(fp, "%s \tserver-url-file = %s\n", prefix, args_info->url_file.c_str());
    fprintf(fp, "%s \tdaemon          = %s\n", prefix, (args_info->daemon_flag)?"true":"false");
    fprintf(fp, "%s \tverbose         = %d\n", prefix, args_info->verbose);
    fprintf(fp, "%s \tlogfile         = %s\n", prefix, args_info->logfile.c_str());

    fprintf(fp, "%s -----------------------------------\n", prefix);
}

static void generate_contact_info(char *myid)
{
    NNTI_url *all_urls=NULL;
    int rank, np;
    char contact_path[1024];
    log_level debug_level = netcdf_debug_level;
    debug_level = LOG_ALL;

    log_debug(netcdf_debug_level, "enter");

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    log_debug(debug_level, "rank (%d)", rank);

    if (rank==0) {
        MPI_Comm_size(MPI_COMM_WORLD, &np);
        all_urls=(NNTI_url *)malloc(np*sizeof(NNTI_url));
    }
    MPI_Gather(myid, sizeof(NNTI_url), MPI_BYTE,
               all_urls, sizeof(NNTI_url), MPI_BYTE,
               0, MPI_COMM_WORLD);
    if (rank==0) {
        char *contact_file=getenv("NETCDF_CONTACT_INFO");
        if (contact_file==NULL) {
            log_error(debug_level, "NETCDF_CONTACT_INFO env var is undefined.");
            free(all_urls);
            return;
        }
//        sprintf(contact_path, "%s.%04d", contact_file, rank);
        sprintf(contact_path, "%s.tmp", contact_file);
        log_debug(debug_level, "creating contact file (%s)", contact_path);
        FILE *f=fopen(contact_path, "w");
        if (f==NULL) {
            perror("fopen");
        }
        for (int i=0;i<np;i++) {
            fprintf(f, "%s\n",
                    all_urls[i]);
        }
//        fprintf(f, "%u@%u@%s@%u\n",
//                myid->nid, myid->pid,
//                myid->hostname, (unsigned int)ntohs(myid->port));
        fclose(f);
        rename(contact_path, contact_file);
        free(all_urls);
    }
    log_debug(netcdf_debug_level, "exit");
}


/**
 * @brief The NSSI xfer-server.
 */
int main(int argc, char **argv)
{
    int rc = NSSI_OK;

    nssi_service netcdf_svc;
    struct netcdf_args args_info;
    log_level debug_level;
    char logfile[1024];
    int rank, np;
    char my_url[NSSI_URL_LEN];

    bool success = true;

    // Initialize arguments
    args_info.daemon_flag = false;

    args_info.verbose = LOG_INFO;
    args_info.logfile = "";

    args_info.server_url = "";
    args_info.url_file = "";

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    /**
     * We make extensive use of the \ref Teuchos::CommandLineProcessor for command-line
     * options to control the behavior of the test code.   To evaluate performance,
     * the "num-trials", "num-reqs", and "len" options control the amount of data transferred
     * between client and server.  The "io-method" selects the type of data transfer.  The
     * server-url specifies the URL of the server.  If running as a server, the server-url
     * provides a recommended URL when initializing the network transport.
     */
    try {

        //out << Teuchos::Teuchos_Version() << std::endl << std::endl;

        // Creating an empty command line processor looks like:
        Teuchos::CommandLineProcessor parser;
        parser.setDocString(
                "This example program demonstrates a simple data-transfer service "
                "built using the NEtwork Scalable Service Interface (Nessie)."
        );

        /* To set and option, it must be given a name and default value.  Additionally,
           each option can be given a help std::string.  Although it is not necessary, a help
           std::string aids a users comprehension of the acceptable command line arguments.
           Some examples of setting command line options are:
         */

        parser.setOption("daemon", "no-daemon", &args_info.daemon_flag, "Run server as a daemon process in the background");

        parser.setOption("verbose",(int*)(&args_info.verbose), "Debug level");
        parser.setOption("logfile", &args_info.logfile, "log file");

        parser.setOption("server-url", &args_info.server_url, "URL client uses to find the server");
        parser.setOption("server-url-file", &args_info.url_file, "File that has URL client uses to find server");


        /* There are also two methods that control the behavior of the
           command line processor.  First, for the command line processor to
           allow an unrecognized a command line option to be ignored (and
           only have a warning printed), use:
         */
        parser.recogniseAllOptions(true);

        /* Second, by default, if the parser finds a command line option it
           doesn't recognize or finds the --help option, it will throw an
           std::exception.  If you want prevent a command line processor from
           throwing an std::exception (which is important in this program since
           we don't have an try/catch around this) when it encounters a
           unrecognized option or help is printed, use:
         */
        parser.throwExceptions(false);

        /* We now parse the command line where argc and argv are passed to
           the parse method.  Note that since we have turned off std::exception
           throwing above we had better grab the return argument so that
           we can see what happened and act accordingly.
         */
        Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn= parser.parse( argc, argv );

        if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {
            return 0;
        }

        if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
            return 1; // Error!

        }

        // Here is where you would use these command line arguments but for this example program
        // we will just print the help message with the new values of the command-line arguments.
        //if (rank == 0)
        //    out << "\nPrinting help message with new values of command-line arguments ...\n\n";

        //parser.printHelpMessage(argv[0],out);

    }

    TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

    log_debug(args_info.verbose, "%d: Finished processing arguments", rank);

    /* initialize and enable logging */
    if (args_info.logfile.empty()) {
        logger_init((log_level)args_info.verbose, NULL);
    } else {
        sprintf(logfile, "%s.%04d", args_info.logfile.c_str(), rank);
        logger_init((log_level)args_info.verbose, logfile);
    }
    netcdf_debug_level=(log_level)(log_level)args_info.verbose;
    debug_level = (log_level)args_info.verbose;

    if (args_info.daemon_flag) {
        nssi_daemonize();
    }

    memset(&netcdf_svc, 0, sizeof(nssi_service));

    nssi_rpc_init(NSSI_DEFAULT_TRANSPORT, NSSI_DEFAULT_ENCODE, NULL);

    nssi_get_url(NSSI_DEFAULT_TRANSPORT, my_url, NSSI_URL_LEN);
    generate_contact_info(my_url);

    log_debug(debug_level, "Initialize netcdf service");

    /* initialize the nssi service */
    rc = nssi_service_init(NSSI_DEFAULT_TRANSPORT, NSSI_SHORT_REQUEST_SIZE, &netcdf_svc);
    if (rc != NSSI_OK) {
        log_error(debug_level, "could not init xfer_svc: %s",
                nssi_err_str(rc));
        return -1;
    }

    /* initialize netcdf service */
    rc = netcdf_server_init();

    /* print the arguments */
    print_opts(stdout, &args_info, "%");

    //mlockall(MCL_FUTURE);

    /* start processing requests */
    netcdf_svc.max_reqs = -1;
    rc = nssi_service_start(&netcdf_svc);
    if (rc != NSSI_OK) {
        log_info(netcdf_debug_level, "exited xfer_svc: %s",
                nssi_err_str(rc));
    }

    /* shutdown the xfer_svc */
    log_debug(debug_level, "shutting down service library");
    nssi_service_fini(&netcdf_svc);

    nssi_rpc_fini(NSSI_DEFAULT_TRANSPORT);

    MPI_Finalize();

    return rc;
}

#endif // HAVE_TRIOS_PNETCDF
