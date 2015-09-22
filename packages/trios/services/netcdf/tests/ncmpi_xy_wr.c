/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This is a very simple example which writes a 2D array of
   sample data. To handle this in netCDF we create two shared
   dimensions, "x" and "y", and a netCDF variable, called "data".

   This example demonstrates the netCDF C API. This is part of the
   netCDF tutorial, which can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

   Full documentation of the netCDF C API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c

   $Id: simple_xy_wr.c,v 1.11 2006/09/27 13:44:36 ed Exp $
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pnetcdf.h>
#include <mpi.h>

/* This is the name of the data file we will create. */
#define FILE_NAME "ncmpi_xy.nc"

/* We are writing 2D data, a 6 x 12 grid. */
#define NDIMS 2
#define NX 6
#define NY 12

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", ncmpi_strerror(e)); exit(ERRCODE);}


int
main(int argc, char *argv[])
{
    /* This will be the netCDF ID for the file and data variable. */
    int ncid, varid;

    int *data_out;
    int data_len[2] = {NX, NY};
    int dimids[2];
    int x_dimid;
    int y_dimid;
    int bufcount;

    /* Loop indexes, and error handling. */
    int i, x, y, retval;
    int proc_dims[2];
    int proc_periods[2];
    int proc_coord[2];
    int proc_ndims = 2;   /* use 2-d topology for processors */

    /* domain variables */
    MPI_Offset start[2];
    MPI_Offset count[2];
    MPI_Offset stride[2];

    /* MPI variables and communicators */
    int np, rank;
    MPI_Comm comm_cart;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    memset(proc_dims, 0, 2*sizeof(int));
    memset(proc_periods, 0, 2*sizeof(int));

    memset(start, 0, 2*sizeof(MPI_Offset));
    memset(count, 0, 2*sizeof(MPI_Offset));
    memset(stride, 0, 2*sizeof(MPI_Offset));

    /* create the dimensions */
    MPI_Dims_create(np, proc_ndims, proc_dims);


    if (rank == 0) {
        fprintf(stdout, "Global Parameters: np=%d, proc_dims = [%d, %d], data_len=[%d, %d]\n\n",
                np, proc_dims[0], proc_dims[1], NX, NY);
    }

    /* create the cartiesion communicator */
    MPI_Cart_create(MPI_COMM_WORLD, proc_ndims, proc_dims, proc_periods, 0, &comm_cart);

    /* get the coordinates of this process */
    MPI_Cart_coords(comm_cart, rank, 2, proc_coord);

    /* Do the domain decomposition ... assign start, offset, count, stide */

    /* X dimension */
    for (i=0; i<proc_ndims; i++) {
        /* calculate number of items for this processor in dimension i */
        count[i] = data_len[i] / proc_dims[i];
        if (proc_coord[i] < (data_len[i] % proc_dims[i]))
            count[i]++;

        /* calculate starting index for this processor in dimension i */
        start[i] = proc_coord[i] * (data_len[i] / proc_dims[i]);
        if (proc_coord[i] < (data_len[i] % proc_dims[i]))
            start[i] += (proc_coord[i]);  /* add 1 for each large node prior to this node */
        else
            start[i] += data_len[i] % proc_dims[i]; /* add one for each large node */

        /* calculate stride for this processor in dimension i */
        stride[i] = data_len[i] - count[i];
        stride[i] = 1;
    }

    /* print domain information in order */
    for (i=0; i<np; i++) {
        if (rank == i) {
            fprintf(stdout, "%d: Domain: "
                    "Coords = [%d, %d], "
                    "Start  = [%d, %d], "
                    "Count  = [%d, %d], "
                    "Stride = [%d, %d]\n", rank,
                    proc_coord[0], proc_coord[1],
                    (int)start[0], (int)start[1], (int)count[0],
                    (int)count[1], (int)stride[0], (int)stride[1]);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* allocate space for data */
    bufcount = count[0]*count[1];
    data_out = (int *)calloc(bufcount, sizeof(int));


    /* Create the file. The NC_CLOBBER parameter tells netCDF to
     * overwrite this file, if it already exists.*/
    if ((retval = ncmpi_create(MPI_COMM_WORLD, FILE_NAME, NC_CLOBBER, MPI_INFO_NULL, &ncid)))
        ERR(retval);

    /* Define the dimensions. NetCDF will hand back an ID for each. */
    if ((retval = ncmpi_def_dim(ncid, "x", NX, &x_dimid)))
        ERR(retval);
    if ((retval = ncmpi_def_dim(ncid, "y", NY, &y_dimid)))
        ERR(retval);

    /* The dimids array is used to pass the IDs of the
     * dimensions of the variable. */
    dimids[0] = x_dimid;
    dimids[1] = y_dimid;

    /* Define the variable. The type of the variable in this case is
     * NC_INT (4-byte integer). */
    if ((retval = ncmpi_def_var(ncid, "data", NC_INT, NDIMS,
            dimids, &varid)))
        ERR(retval);

    /* End define mode. This tells netCDF we are done defining
     * metadata. */
    if ((retval = ncmpi_enddef(ncid)))
        ERR(retval);


    /* set the data. */
    for (x = 0; x < count[0]; x++) {
        for (y = 0; y < count[1]; y++) {
            int index = x*count[1] + y;
            int realx = start[0] + x;
            int realy = start[1] + y;
            data_out[index] = realx * NY + realy;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);


    /* Put the application in independent mode */
    if ((retval = ncmpi_begin_indep_data(ncid)))
        ERR(retval);

    /* Read the data. */
    if ((retval = ncmpi_put_vars(ncid, varid, start, count, stride, data_out, bufcount, MPI_INT))) {

        fprintf(stderr, "%d: %s\n", rank, ncmpi_strerror(retval));
        goto cleanup;
    }

    /* Put the application in independent mode */
    if ((retval = ncmpi_end_indep_data(ncid)))
        ERR(retval);

    /* Close the file, freeing all resources. */
    if ((retval = ncmpi_close(ncid)))
        ERR(retval);

    if (rank == 0)
        printf("\n*** SUCCESS writing example file %s!\n", FILE_NAME);


cleanup:
    MPI_Finalize();
    if (data_out) free(data_out);

    return 0;
}
