/***********************************************************
 *
 * This test program reads a netCDF file, and then write it
 * out to another netCDF file, using the parallel netCDF
 * library using MPI-IO. The two files should be identical.
 *
 * Input File: "../data/test_int.nc"  generated from original netcdf-3.
 * Output File: "testread.nc"
 *
 *  The CDL notation for the test dataset is shown below:
 *
 *    netcdf test {
 *
 *       dimensions:
 *
 *            x = 100, y = 100, z = 100, time = NC_UNLIMITED;
 *
 *
 *       variables:  // variable types, names, shapes, attributes
 *
 *            int   square(x, y);
 *                     squre: description = "2-D integer array";
 *
 *            int   cube(x,y,z);
 *
 *            int   time(time);  // coordinate & record variable
 *
 *            int   xytime(time, x, y);  // record variable
 *
 *
 *      // global attributes
 *
 *           :title = "example netCDF dataset";
 *
 *
 *      data:  // data written for variables
 *          square  = 0, 1, 2, 3,  ... , 9999;
 *          cube    = 0, 1, 2, 3,  ... , 999999;
 *          time    = 0, 1, 2, 3,  ... , 99;    // 100 records
 *          xytime  = 0, 1, 2, 3,  ... , 9999;  // 100 records
 *   }
 *
 *
 *
 * This test uses collective APIs to read/write variable data and
 * only deals with integer variables.
 *
 * This test assume # of processors = 4
 *
 **********************************************************/

#include <stdio.h>
#include <mpi.h>
#include <pnetcdf.h>
#include <stdlib.h>
#include "testutils.h"

/* Prototype for functions used only in this file */
static void handle_error(int status);

static void handle_error(int status) {
  printf("%s\n", ncmpi_strerror(status));
}

int main(int argc, char **argv) {

  int i, j;
  int status;
  int ncid1, ncid2;
  int ndims, nvars, ngatts, unlimdimid;
  char name[NC_MAX_NAME];
  nc_type type, vartypes[NC_MAX_VARS];
  MPI_Offset attlen;
  MPI_Offset dimlen, shape[NC_MAX_VAR_DIMS], varsize, start[NC_MAX_VAR_DIMS];
  void *valuep;
  int dimids[NC_MAX_DIMS], varids[NC_MAX_VARS];
  int vardims[NC_MAX_VARS][NC_MAX_VAR_DIMS/16]; /* divided by 16 due to my memory limitation */
  int varndims[NC_MAX_VARS], varnatts[NC_MAX_VARS];
  params opts;

  int rank;
  int nprocs;
  MPI_Comm comm = MPI_COMM_WORLD;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
          fprintf(stderr, "Testing read ... ");
  parse_read_args(argc, argv, rank, &opts);

  /**********  START OF NETCDF ACCESS **************/


  /* Read a netCDF file and write it out to another file */

  /**
   * Open the input dataset - ncid1:
   *   File name: "../data/test_int.nc"
   *   Dataset API: Collective
   * And create the output dataset - ncid2:
   *   File name: "testread.nc"
   *   Dataset API: Collective
   */

  status = ncmpi_open(comm, opts.infname, 0, MPI_INFO_NULL, &ncid1);
  if (status != NC_NOERR) handle_error(status);

  status = ncmpi_create(comm, opts.outfname, NC_CLOBBER, MPI_INFO_NULL, &ncid2);
  if (status != NC_NOERR) handle_error(status);


  /**
   * Inquire the dataset definitions of input dataset AND
   * Add dataset definitions for output dataset.
   */

  status = ncmpi_inq(ncid1, &ndims, &nvars, &ngatts, &unlimdimid);
  if (status != NC_NOERR) handle_error(status);


  /* Inquire global attributes, assume CHAR attributes. */

  for (i = 0; i < ngatts; i++) {
    status = ncmpi_inq_attname(ncid1, NC_GLOBAL, i, name);
    if (status != NC_NOERR) handle_error(status);
    status = ncmpi_inq_att (ncid1, NC_GLOBAL, name, &type, &attlen);
    if (status != NC_NOERR) handle_error(status);
    switch (type) {
      case NC_CHAR:
        valuep = (void *)malloc(attlen * sizeof(char));
        status = ncmpi_get_att_text(ncid1, NC_GLOBAL, name, valuep);
        if (status != NC_NOERR) handle_error(status);
        status = ncmpi_put_att_text (ncid2, NC_GLOBAL, name, attlen, (char *)valuep);
        if (status != NC_NOERR) handle_error(status);
        free(valuep);
        break;
      case NC_SHORT:
        valuep = (void *)malloc(attlen * sizeof(short));
        status = ncmpi_get_att_short(ncid1, NC_GLOBAL, name, valuep);
        if (status != NC_NOERR) handle_error(status);
        status = ncmpi_put_att_short (ncid2, NC_GLOBAL, name, type, attlen, (short *)valuep);
        if (status != NC_NOERR) handle_error(status);
        free(valuep);
        break;
      case NC_INT:
        valuep = (void *)malloc(attlen * sizeof(int));
        status = ncmpi_get_att_int(ncid1, NC_GLOBAL, name, valuep);
        if (status != NC_NOERR) handle_error(status);
        status = ncmpi_put_att_int (ncid2, NC_GLOBAL, name, type, attlen, (int *)valuep);
        if (status != NC_NOERR) handle_error(status);
        free(valuep);
        break;
      case NC_FLOAT:
        valuep = (void *)malloc(attlen * sizeof(float));
        status = ncmpi_get_att_float(ncid1, NC_GLOBAL, name, valuep);
        if (status != NC_NOERR) handle_error(status);
        status = ncmpi_put_att_float (ncid2, NC_GLOBAL, name, type, attlen, (float *)valuep);
        if (status != NC_NOERR) handle_error(status);
        free(valuep);
        break;
      case NC_DOUBLE:
        valuep = (void *)malloc(attlen * sizeof(double));
        status = ncmpi_get_att_double(ncid1, NC_GLOBAL, name, valuep);
        if (status != NC_NOERR) handle_error(status);
        status = ncmpi_put_att_double (ncid2, NC_GLOBAL, name, type, attlen, (double *)valuep);
        if (status != NC_NOERR) handle_error(status);
        free(valuep);
        break;
      default:
        ;
        /* handle unexpected types */
    }
  }

  /* Inquire dimension */

  for (i = 0; i < ndims; i++) {
    status = ncmpi_inq_dim(ncid1, i, name, &dimlen);
    if (status != NC_NOERR) handle_error(status);
    if (i == unlimdimid)
      dimlen = NC_UNLIMITED;
    status = ncmpi_def_dim(ncid2, name, dimlen, dimids+i);
    if (status != NC_NOERR) handle_error(status);
  }

  /* Inquire variables */

  for (i = 0; i < nvars; i++) {
    status = ncmpi_inq_var (ncid1, i, name, vartypes+i, varndims+i, vardims[i], varnatts+i);
    if (status != NC_NOERR) handle_error(status);

    status = ncmpi_def_var(ncid2, name, vartypes[i], varndims[i], vardims[i], varids+i);
    if (status != NC_NOERR) handle_error(status);

    /* var attributes, assume CHAR attributes */

    for (j = 0; j < varnatts[i]; j++) {
      status = ncmpi_inq_attname(ncid1, varids[i], j, name);
      if (status != NC_NOERR) handle_error(status);
      status = ncmpi_inq_att (ncid1, varids[i], name, &type, &attlen);
      if (status != NC_NOERR) handle_error(status);
      switch (type) {
        case NC_CHAR:
          valuep = (void *)malloc(attlen * sizeof(char));
          status = ncmpi_get_att_text(ncid1, varids[i], name, valuep);
          if (status != NC_NOERR) handle_error(status);
          status = ncmpi_put_att_text (ncid2, varids[i], name, attlen, (char *)valuep);
          if (status != NC_NOERR) handle_error(status);
          free(valuep);
          break;
        case NC_SHORT:
          valuep = (void *)malloc(attlen * sizeof(short));
          status = ncmpi_get_att_short(ncid1, varids[i], name, valuep);
          if (status != NC_NOERR) handle_error(status);
          status = ncmpi_put_att_short (ncid2, varids[i], name, type, attlen, (short *)valuep);
          if (status != NC_NOERR) handle_error(status);
          free(valuep);
          break;
        case NC_INT:
          valuep = (void *)malloc(attlen * sizeof(int));
          status = ncmpi_get_att_int(ncid1, varids[i], name, valuep);
          if (status != NC_NOERR) handle_error(status);
          status = ncmpi_put_att_int (ncid2, varids[i], name, type, attlen, (int *)valuep);
          if (status != NC_NOERR) handle_error(status);
          free(valuep);
          break;
        case NC_FLOAT:
          valuep = (void *)malloc(attlen * sizeof(float));
          status = ncmpi_get_att_float(ncid1, varids[i], name, valuep);
          if (status != NC_NOERR) handle_error(status);
          status = ncmpi_put_att_float (ncid2, varids[i], name, type, attlen, (float *)valuep);
          if (status != NC_NOERR) handle_error(status);
          free(valuep);
          break;
        case NC_DOUBLE:
          valuep = (void *)malloc(attlen * sizeof(double));
          status = ncmpi_get_att_double(ncid1, varids[i], name, valuep);
          if (status != NC_NOERR) handle_error(status);
          status = ncmpi_put_att_double (ncid2, varids[i], name, type, attlen, (double *)valuep);
          if (status != NC_NOERR) handle_error(status);
          free(valuep);
          break;
        default:
          ;
        /* handle unexpected types */
      }
    }
  }

  /**
   * End Define Mode (switch to data mode) for output dataset
   *   Dataset API: Collective
   */

  status = ncmpi_enddef(ncid2);
  if (status != NC_NOERR) handle_error(status);

  /**
   * Read data of variables from input dataset
   * (ONLY DEAL WITH: NC_INT, NC_FLOAT, NC_DOUBLE for now)
   * Write the data out to the corresponding variables in the output dataset
   *
   *  Data Partition (Assume 4 processors):
   *   square: 2-D, (Block, *), 25*100 from 100*100
   *   cube:   3-D, (Block, *, *), 25*100*100 from 100*100*100
   *   xytime: 3-D, (Block, *, *), 25*100*100 from 100*100*100
   *   time:   1-D, Block-wise, 25 from 100
   *
   *  Data Mode API: collective
   */

  for (i = 0; i < NC_MAX_VAR_DIMS; i++)
    start[i] = 0;
  for (i = 0; i < nvars; i++) {
    varsize = 1;
    for (j = 0; j < varndims[i]; j++) {
      status = ncmpi_inq_dim(ncid1, vardims[i][j], name, shape + j);
      if (status != NC_NOERR) handle_error(status);
      if (j == 0) {
        shape[j] /= nprocs;
        start[j] = shape[j] * rank;
      }
      varsize *= shape[j];
    }
    switch (vartypes[i]) {
      case NC_CHAR:
        break;
      case NC_SHORT:
        valuep = (void *)malloc(varsize * sizeof(short));
        status = ncmpi_get_vara_short_all(ncid1, i, start, shape, (short *)valuep);
        if (status != NC_NOERR) handle_error(status);
        status = ncmpi_put_vara_short_all(ncid2, varids[i],
                                     start, shape, (short *)valuep);
        if (status != NC_NOERR) handle_error(status);
        free(valuep);
        break;
      case NC_INT:
        valuep = (void *)malloc(varsize * sizeof(int));
        status = ncmpi_get_vara_int_all(ncid1, i, start, shape, (int *)valuep);
        if (status != NC_NOERR) handle_error(status);
        status = ncmpi_put_vara_int_all(ncid2, varids[i],
                                     start, shape, (int *)valuep);
        if (status != NC_NOERR) handle_error(status);
        free(valuep);
        break;
      case NC_FLOAT:
        valuep = (void *)malloc(varsize * sizeof(float));
        status = ncmpi_get_vara_float_all(ncid1, i, start, shape, (float *)valuep);
        if (status != NC_NOERR) handle_error(status);
        status = ncmpi_put_vara_float_all(ncid2, varids[i],
                                     start, shape, (float *)valuep);
        if (status != NC_NOERR) handle_error(status);
        free(valuep);
        break;
      case NC_DOUBLE:
        valuep = (void *)malloc(varsize * sizeof(double));
        status = ncmpi_get_vara_double_all(ncid1, i, start, shape, (double *)valuep);
        if (status != NC_NOERR) handle_error(status);
        status = ncmpi_put_vara_double_all(ncid2, varids[i],
                                     start, shape, (double *)valuep);
        if (status != NC_NOERR) handle_error(status);
        free(valuep);
        break;
      default:
        ;
        /* handle unexpected types */
    }
  }

  status = ncmpi_sync(ncid1);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_sync(ncid2);
  if (status != NC_NOERR) handle_error(status);

  /**
   * Close the datasets
   *   Dataset API:  collective
   */

  status = ncmpi_close(ncid1);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_close(ncid2);
  if (status != NC_NOERR) handle_error(status);

  /*******************  END OF NETCDF ACCESS  ****************/

if (rank == 0)
  fprintf(stderr, "OK\nInput file %s copied to: %s!\n", opts.infname, opts.outfname);

  MPI_Finalize();
  return 0;
}
