/* This is part of the netCDF package.
   Copyright 2006 University Corporation for Atmospheric Research/Unidata.
   See COPYRIGHT file for conditions of use.

   This is a simple example which reads a small dummy array, which was
   written by simple_xy_wr.c. This is intended to illustrate the use
   of the netCDF C API.

   This program is part of the netCDF tutorial:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

   Full documentation of the netCDF C API can be found at:
   http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c

   $Id: simple_xy_rd.c,v 1.9 2006/08/17 23:00:55 russ Exp $
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pnetcdf.h>
#include <mpi.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "ncmpi_xy.nc"

/* We are reading 2D data, a 6 x 12 grid. */
#define NX 6
#define NY 12

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", ncmpi_strerror(e)); goto cleanup;}

int
main(int argc, char *argv[])
{
   /* This will be the netCDF ID for the file and data variable. */
   int ncid, varid;

   int *data_in = NULL;
   int data_len[2] = {NX, NY};

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


   /* Open the file. NC_NOWRITE tells netCDF we want read-only access
    * to the file.*/
   if ((retval = ncmpi_open(MPI_COMM_WORLD, FILE_NAME, NC_NOWRITE, MPI_INFO_NULL, &ncid)))
      ERR(retval);

   /* Get the varid of the data variable, based on its name. */
   if ((retval = ncmpi_inq_varid(ncid, "data", &varid)))
      ERR(retval);


   /* Put the application in independent mode */
   if ((retval = ncmpi_begin_indep_data(ncid)))
         ERR(retval);

   /* Read the data. */
   int bufcount = count[0] * count[1];
   data_in = (int *)calloc(bufcount, sizeof(int));
   if ((retval = ncmpi_get_vars(ncid, varid, start, count, stride, data_in, bufcount, MPI_INT))) {

       fprintf(stderr, "%d: %s\n", rank, ncmpi_strerror(retval));
       goto cleanup;
   }


   /* Put the application in independent mode */
   if ((retval = ncmpi_end_indep_data(ncid)))
       ERR(retval);

   /* set the data. */
   for (x = 0; x < count[0]; x++) {
       for (y = 0; y < count[1]; y++) {
           int index = x*count[1] + y;
           int realx = start[0] + x;
           int realy = start[1] + y;
           data_in[index] = realx * NY + realy;
       }
   }


   /* Close the file, freeing all resources. */
   if ((retval = ncmpi_close(ncid)))
      ERR(retval);

   if (rank == 0)
       printf("\n*** SUCCESS reading example file %s!\n", FILE_NAME);


cleanup:
   MPI_Finalize();
   if (data_in) free(data_in);

   return 0;
}
