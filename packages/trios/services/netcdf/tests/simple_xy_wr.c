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
#include <netcdf.h>

/* This is the name of the data file we will create. */
#define FILE_NAME "simple_xy.nc"

/* We are writing 2D data, a 6 x 12 grid. */
#define NDIMS 2
#define NX 6
#define NY 12

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int
main()
{
   /* When we create netCDF variables and dimensions, we get back an
    * ID for each one. */
   int ncid, x_dimid, y_dimid, varid;
   int dimids[NDIMS];

   /* This is the data array we will write. It will just be filled
    * with a progression of numbers for this example. */
   int data_out[NX][NY];

   /* Loop indexes, and error handling. */
   int x, y, retval;

   /* Create some pretend data. If this wasn't an example program, we
    * would have some real data to write, for example, model
    * output. */
   for (x = 0; x < NX; x++)
      for (y = 0; y < NY; y++)
	 data_out[x][y] = x * NY + y;

   /* Always check the return code of every netCDF function call. In
    * this example program, any retval which is not equal to NC_NOERR
    * (0) will cause the program to print an error message and exit
    * with a non-zero return code. */

   /* Create the file. The NC_CLOBBER parameter tells netCDF to
    * overwrite this file, if it already exists.*/
   if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &ncid)))
      ERR(retval);

   /* Define the dimensions. NetCDF will hand back an ID for each. */
   if ((retval = nc_def_dim(ncid, "x", NX, &x_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid, "y", NY, &y_dimid)))
      ERR(retval);

   /* The dimids array is used to pass the IDs of the dimensions of
    * the variable. */
   dimids[0] = x_dimid;
   dimids[1] = y_dimid;

   /* Define the variable. The type of the variable in this case is
    * NC_INT (4-byte integer). */
   if ((retval = nc_def_var(ncid, "data", NC_INT, NDIMS, 
			    dimids, &varid)))
      ERR(retval);

   /* End define mode. This tells netCDF we are done defining
    * metadata. */
   if ((retval = nc_enddef(ncid)))
      ERR(retval);

   /* Write the pretend data to the file. Although netCDF supports
    * reading and writing subsets of data, in this case we write all
    * the data in one operation. */
   if ((retval = nc_put_var_int(ncid, varid, &data_out[0][0])))
      ERR(retval);

   /* Close the file. This frees up any internal netCDF resources
    * associated with the file, and flushes any buffers. */
   if ((retval = nc_close(ncid)))
      ERR(retval);

   printf("*** SUCCESS writing example file simple_xy.nc!\n");
   return 0;
}
