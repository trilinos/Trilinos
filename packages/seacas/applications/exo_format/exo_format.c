/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/*
 * Determine the format of the exodus file:
 * netcdf format: classic, 64-bit offset, netcdf-4, netcdf-4-classic
 * exodus version,
 * integer sizes for ids, maps, bulk data
 * exodus file-size attribute (normal, large)
 */

#include "exodusII.h"
#include "netcdf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#define NOT_NETCDF 2
#define NOT_EXODUSII 3

char *progname;
char *filename;

int main(int argc, char *argv[])
{
  int   c1;
  int   c2;
  int   c3;
  int   c4;
  FILE *fid;

  size_t j;
  size_t k;

  int   exoid;
  int   CPU_word_size;
  int   IO_word_size;
  float version;
  int   file_size;
  int   netcdf_based    = 0;
  int   hdf5_based      = 0;
  int   nc_format       = 0;
  int   int64_status    = 0;
  int   max_name_length = 0;
  int   fn_idx          = 1;
  char  cversion[9];

  CPU_word_size = 0; /* float or double */
  IO_word_size  = 0; /* use what is stored in file */

  /* Determine if filename was given */
  progname = argv[0];

  if (argc <= 1) {
    fprintf(stderr, "USAGE: %s [-config] {filename}\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (argv[1][0] == '-') {
    if (argv[1][1] == 'c') {
      fprintf(stderr, "\nExodus Configuration Information:\n");
      ex_print_config();
    }
    fn_idx = 2;
    if (argc <= 2) {
      exit(0);
    }
  }

  /* examine file */
  filename = argv[fn_idx]; /* filename path */

  fid = fopen(filename, "r");
  if (fid == NULL) {
    (void)fprintf(stderr, "         Could not open %s\n", filename);
    exit(EXIT_FAILURE);
  }
  c1 = getc(fid);
  c2 = getc(fid);
  c3 = getc(fid);
  c4 = getc(fid);
  fclose(fid);
  if (c1 == 'C' && c2 == 'D' && c3 == 'F') {
    netcdf_based = 1;
  }
  else if (c2 == 'H' && c3 == 'D' && c4 == 'F') {
    hdf5_based = 1;
  }
  else {
    (void)fprintf(stderr, "         %s is not an EXODUS or netCDF file\n", filename);
    exit(NOT_NETCDF);
  }
  exoid = ex_open(filename, EX_READ, /* access mode = READ */
                  &CPU_word_size,    /* CPU word size */
                  &IO_word_size,     /* IO word size */
                  &version);         /* Exodus library version */

  if (exoid < 0) {
    if (netcdf_based) {
      (void)fprintf(stderr, "         %s is a NetCDF file, but not a valid EXODUS file\n",
                    filename);
    }
    else if (hdf5_based) {
      (void)fprintf(stderr, "         %s is an HDF5 file, but not a valid EXODUS file.\n",
                    filename);
    }
    exit(NOT_EXODUSII);
  }

  file_size = ex_large_model(exoid);

  fprintf(stderr, "\n\t%s is an EXODUS file, version %4.2f\n\n", filename, version);

  /* Determine int sizes */
  int64_status = ex_int64_status(exoid);
  if (int64_status & EX_IDS_INT64_DB) {
    fprintf(stderr, "\t\tIDs are stored as 64-bit integers\n");
  }
  else {
    fprintf(stderr, "\t\tIDs are stored as 32-bit integers\n");
  }

  if (int64_status & EX_MAPS_INT64_DB) {
    fprintf(stderr, "\t\tMap entries are stored as 64-bit integers\n");
  }
  else {
    fprintf(stderr, "\t\tMap entries are stored as 32-bit integers\n");
  }

  if (int64_status & EX_BULK_INT64_DB) {
    fprintf(stderr, "\t\tBulk data are stored as 64-bit integers\n");
  }
  else {
    fprintf(stderr, "\t\tBulk data are stored as 32-bit integers\n");
  }

  if (IO_word_size == 4) {
    fprintf(stderr, "\t\tFloating point data are stored as 32-bit floats\n");
  }
  else {
    fprintf(stderr, "\t\tFloating point data are stored as 64-bit doubles\n");
  }
  max_name_length = (int)ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  fprintf(stderr, "\n\t\tMaximum name length is %d\n\n", max_name_length);

  if (file_size == 0) {
    fprintf(stderr, "\t\tFile size attribute is 'normal model'\n");
  }
  else {
    fprintf(stderr, "\t\tFile size attribute is 'large model'\n");
  }

  /* Determine netcdf file version... */
  nc_inq_format(exoid, &nc_format);
  if (nc_format == NC_FORMAT_CLASSIC) {
    fprintf(stderr, "\t\tNetCDF Variant is 'classic'\n");
  }
  else if (nc_format == NC_FORMAT_64BIT) {
    fprintf(stderr, "\t\tNetCDF Variant is '64-bit offset'\n");
  }
#if defined NC_FORMAT_64BIT_DATA
  else if (nc_format == NC_FORMAT_64BIT_DATA) {
    fprintf(stderr, "\t\tNetCDF Variant is '64-bit data (CDF5)'\n");
  }
#endif
  else if (nc_format == NC_FORMAT_NETCDF4) {
    fprintf(stderr, "\t\tNetCDF Variant is 'netCDF-4'\n");
  }
  else if (nc_format == NC_FORMAT_NETCDF4_CLASSIC) {
    fprintf(stderr, "\t\tNetCDF Variant is 'netCDF-4 classic'\n");
  }
  else {
    fprintf(stderr, "\t\tNetCDF Variant is 'unknown'\n");
  }

  if (netcdf_based) {
    fprintf(stderr, "\t\tUnderlying data format is netcdf\n");
  }
  if (hdf5_based) {
    fprintf(stderr, "\t\tUnderlying data format is hdf5\n");
  }

  fprintf(stderr, "\n");

  /* Determine number of dims and vars -- useful in debugging incorrect NC_MAX_DIMS|VARS in netcdf.h
   */
#if NC_HAS_NC4
  {
    int ndims = 0;
    int nvars = 0;
    nc_inq_dimids(exoid, &ndims, NULL, 0);
    nc_inq_varids(exoid, &nvars, NULL);
    fprintf(stderr, "\t\tNumber of dims = %d\n", ndims);
    fprintf(stderr, "\t\tNumber of vars = %d\n", nvars);
  }
#endif

  if (ex_close(exoid) == -1) {
    printf("ex_close failed");
  }

  version += 0.00005F;
  sprintf(cversion, "%4.2f", version);

  k = strlen(cversion);
  for (j = 0; j < k; j++) {
    if (cversion[j] == '.') {
      break;
    }
  }
  if (j == k) {
    (void)fprintf(stderr, "         %s is not an EXODUS file\n", filename);
    exit(NOT_EXODUSII);
  }

  exit(-2);
}
