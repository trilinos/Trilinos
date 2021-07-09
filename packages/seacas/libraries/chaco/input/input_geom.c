/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "defs.h"    // for FALSE, TRUE
#include "smalloc.h" // for smalloc
#include <stdio.h>   // for printf, fclose, fscanf, EOF, etc

int input_geom(FILE *  fingeom,                /* geometry input file */
               char *  geomname,               /* name of geometry file */
               int     nvtxs,                  /* number of coordinates to read */
               int *   igeom,                  /* dimensionality of geometry */
               float **x, float **y, float **z /* coordinates of vertices */
)
{
  extern int CHECK_INPUT; /* print any warning messages? */
  extern int DEBUG_INPUT; /* echo that read was successful? */
  extern int DEBUG_TRACE; /* trace main execution path */
  float      xc, yc, zc;  /* first x, y, z coordinate */
  int        nread;       /* number of lines of coordinates read */
  int        flag;        /* any bad data at end of file? */
  int        line_num;    /* counts input lines in file */
  int        end_flag;    /* return conditional */
  int        ndims;       /* number of values in an input line */
  int        i = 0;       /* loop counter */
  double     read_val();

  if (DEBUG_TRACE > 0) {
    printf("<Entering input_geom>\n");
  }

  *x = *y = *z = NULL;
  line_num     = 0;
  end_flag     = 1;
  while (end_flag == 1) {
    xc = read_val(fingeom, &end_flag);
    ++line_num;
  }

  if (end_flag == -1) {
    printf("No values found in geometry file `%s'\n", geomname);
    fclose(fingeom);
    return (1);
  }

  ndims = 1;
  yc    = read_val(fingeom, &end_flag);
  if (end_flag == 0) {
    ndims = 2;
    zc    = read_val(fingeom, &end_flag);
    if (end_flag == 0) {
      ndims = 3;
      read_val(fingeom, &end_flag);
      if (!end_flag) {
        printf("Too many values on input line of geometry file `%s'\n", geomname);

        printf(" Maximum dimensionality is 3\n");
        fclose(fingeom);
        return (1);
      }
    }
  }

  *igeom = ndims;

  *x      = smalloc(nvtxs * sizeof(float));
  (*x)[0] = xc;
  if (ndims > 1) {
    *y      = smalloc(nvtxs * sizeof(float));
    (*y)[0] = yc;
  }
  if (ndims > 2) {
    *z      = smalloc(nvtxs * sizeof(float));
    (*z)[0] = zc;
  }

  for (nread = 1; nread < nvtxs; nread++) {
    ++line_num;
    if (ndims == 1) {
      i = fscanf(fingeom, "%f", &((*x)[nread]));
    }
    else if (ndims == 2) {
      i = fscanf(fingeom, "%f%f", &((*x)[nread]), &((*y)[nread]));
    }
    else if (ndims == 3) {
      i = fscanf(fingeom, "%f%f%f", &((*x)[nread]), &((*y)[nread]), &((*z)[nread]));
    }

    if (i == EOF) {
      printf("Too few lines of values in geometry file; nvtxs=%d, but only %d read\n", nvtxs,
             nread + 1);
      fclose(fingeom);
      return (1);
    }
    if (i != ndims) {
      printf("Wrong number of values in line %d of geometry file `%s'\n", line_num, geomname);
      fclose(fingeom);
      return (1);
    }
  }

  /* Check for spurious extra stuff in file. */
  flag     = FALSE;
  end_flag = 0;
  while (!flag && end_flag != -1) {
    read_val(fingeom, &end_flag);
    if (!end_flag) {
      flag = TRUE;
    }
  }
  if (flag && CHECK_INPUT) {
    printf("Warning: possible error in geometry file `%s'\n", geomname);
    printf(" Numerical data found after expected end of file\n");
  }

  fclose(fingeom);

  if (DEBUG_INPUT > 0) {
    printf("Finished reading geometry file `%s'; dimension = %d.\n", geomname, ndims);
  }
  return (0);
}
