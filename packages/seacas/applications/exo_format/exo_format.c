/*
 * Copyright(C) 1999-2021, 2023, 2025 National Technology & Engineering Solutions
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
 * change-set (netcdf groups) count
 */

#include "exodusII.h"
#include "netcdf.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#define NOT_NETCDF   2
#define NOT_EXODUSII 3

int main(int argc, char *argv[])
{
  /* Determine if filename was given */
  if (argc <= 1) {
    fprintf(stderr, "USAGE: %s [-config] {filename}\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  int fn_idx = 1;
  if (argv[1][0] == '-') {
    if (argv[1][1] == 'c' || (argv[1][1] == '-' && argv[1][2] == 'c')) {
      fprintf(stderr, "\nExodus Configuration Information:\n");
      ex_print_config();
    }
    fn_idx = 2;
    if (argc <= 2) {
      exit(0);
    }
  }

  /* examine file */
  char *filename = argv[fn_idx]; /* filename path */

  FILE *fid = fopen(filename, "r");
  if (fid == NULL) {
    fprintf(stderr, "         Could not open %s\n", filename);
    exit(EXIT_FAILURE);
  }
  int c1 = getc(fid);
  int c2 = getc(fid);
  int c3 = getc(fid);
  int c4 = getc(fid);
  fclose(fid);

  bool netcdf_based = false;
  bool hdf5_based   = false;
  if (c1 == 'C' && c2 == 'D' && c3 == 'F') {
    netcdf_based = true;
  }
  else if (c2 == 'H' && c3 == 'D' && c4 == 'F') {
    hdf5_based = true;
  }
  else {
    fprintf(stderr, "         %s is not an EXODUS or netCDF file\n", filename);
    exit(NOT_NETCDF);
  }
  float version;
  int   CPU_word_size = 0; /* float or double */
  int   IO_word_size  = 0; /* use what is stored in file */

  int exoid = ex_open(filename, EX_READ, /* access mode = READ */
                      &CPU_word_size,    /* CPU word size */
                      &IO_word_size,     /* IO word size */
                      &version);         /* Exodus library version */

  if (exoid < 0) {
    if (netcdf_based) {
      fprintf(stderr, "         %s is a NetCDF file, but not a valid EXODUS file\n", filename);
    }
    else if (hdf5_based) {
      fprintf(stderr, "         %s is an HDF5 file, but not a valid EXODUS file.\n", filename);
    }
    exit(NOT_EXODUSII);
  }

  int file_size = ex_large_model(exoid);

  fprintf(stderr, "\n\t%s is an EXODUS file, version %4.2f\n\n", filename, version);

  /* Determine int sizes */
  int int64_status = ex_int64_status(exoid);
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
  int max_name_length = (int)ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  fprintf(stderr, "\n\t\tMaximum name length is %d\n", max_name_length);

  int num_change_sets = ex_inquire_int(exoid, EX_INQ_NUM_CHILD_GROUPS);
  if (num_change_sets > 0) {
    fprintf(stderr, "\t\tFile contains %d change sets (groups)\n", num_change_sets);
  }

  if (file_size == 0) {
    fprintf(stderr, "\n\t\tFile size attribute is 'normal model'\n");
  }
  else {
    fprintf(stderr, "\n\t\tFile size attribute is 'large model'\n");
  }

  /* Determine netcdf file version... */
  int nc_format = 0;
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
  {
    int ndims = 0;
    int nvars = 0;
#if NC_HAS_HDF5
    nc_inq_dimids(exoid, &ndims, NULL, 0);
    nc_inq_varids(exoid, &nvars, NULL);
#else
    nc_inq(exoid, &ndims, &nvars, NULL, NULL);
#endif
    fprintf(stderr, "\t\tNumber of dims = %d\n", ndims);
    fprintf(stderr, "\t\tNumber of vars = %d\n\n", nvars);
  }

  if (ex_close(exoid) == -1) {
    printf("ex_close failed");
  }

  version += 0.00005F;
  char cversion[9];
  snprintf(cversion, 9, "%4.2f", version);

  size_t k = strlen(cversion);
  size_t j = 0;
  for (; j < k; j++) {
    if (cversion[j] == '.') {
      break;
    }
  }
  if (j == k) {
    fprintf(stderr, "         %s is not an EXODUS file\n", filename);
    exit(NOT_EXODUSII);
  }

  exit(-2);
}
