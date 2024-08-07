// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER



/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ch_input_const.h"
#include "dr_const.h"
#include "dr_externs.h"
#include "dr_compress_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

int chaco_input_geom(
ZOLTAN_FILE *fingeom,		/* geometry input file */
char     *geomname,		/* name of geometry file */
int       nvtxs,		/* number of coordinates to read */
int      *igeom,		/* dimensionality of geometry */
float   **x,         		/* coordinates of vertices */
float   **y,
float   **z
)
{
    const char     *yo = "chaco_input_geom";
    float     xc, yc, zc =0;	/* first x, y, z coordinate */
    int       nread;		/* number of lines of coordinates read */
    int       flag;		/* any bad data at end of file? */
    int       line_num;		/* counts input lines in file */
    int       end_flag;		/* return conditional */
    int       ndims;		/* number of values in an input line */
    int       i=0;		/* loop counter */

    DEBUG_TRACE_START(0, yo);

    *x = *y = *z = NULL;
    line_num = 0;
    end_flag = 1;
    while (end_flag == 1) {
	xc = read_val(fingeom, &end_flag);
	++line_num;
    }

    if (end_flag == -1) {
	printf("No values found in geometry file `%s'\n", geomname);
	ZOLTAN_FILE_close(fingeom);
	return (1);
    }

    ndims = 1;
    yc = read_val(fingeom, &end_flag);
    if (end_flag == 0) {
	ndims = 2;
	zc = read_val(fingeom, &end_flag);
	if (end_flag == 0) {
	    ndims = 3;
	    read_val(fingeom, &end_flag);
	    if (!end_flag) {
		printf("Too many values on input line of geometry file `%s'\n",
		       geomname);

		printf(" Maximum dimensionality is 3\n");
		ZOLTAN_FILE_close(fingeom);
		return (1);
	    }
	}
    }

    *igeom = ndims;

    *x = (float *) malloc((unsigned) nvtxs * sizeof(float));
    (*x)[0] = xc;
    if (ndims > 1) {
	*y = (float *) malloc((unsigned) nvtxs * sizeof(float));
	(*y)[0] = yc;
    }
    if (ndims > 2) {
	*z = (float *) malloc((unsigned) nvtxs * sizeof(float));
	(*z)[0] = zc;
    }

    for (nread = 1; nread < nvtxs; nread++) {
	++line_num;
	if (ndims == 1) {
	    i = ZOLTAN_FILE_scanf(fingeom, "%f", &((*x)[nread]));
	}
	else if (ndims == 2) {
	    i = ZOLTAN_FILE_scanf(fingeom, "%f%f", &((*x)[nread]), &((*y)[nread]));
	}
	else if (ndims == 3) {
	    i = ZOLTAN_FILE_scanf(fingeom, "%f%f%f", &((*x)[nread]), &((*y)[nread]),
		       &((*z)[nread]));
	}

	if (i == EOF) {
	    printf("Too few lines of values in geometry file; nvtxs=%d, but only %d read\n",
		   nvtxs, nread);
	    ZOLTAN_FILE_close(fingeom);
	    return (1);
	}
	else if (i != ndims) {
	    printf("Wrong number of values in line %d of geometry file `%s'\n",
		   line_num, geomname);
	    ZOLTAN_FILE_close(fingeom);
	    return (1);
	}
    }

    /* Check for spurious extra stuff in file. */
    flag = 0;
    end_flag = 0;
    while (!flag && end_flag != -1) {
	read_val(fingeom, &end_flag);
	if (!end_flag)
	    flag = 1;
    }
    if (flag && Debug_Chaco_Input) {
	printf("Warning: possible error in geometry file `%s'\n", geomname);
	printf(" Numerical data found after expected end of file\n");
    }

    ZOLTAN_FILE_close(fingeom);

    DEBUG_TRACE_END(0, yo);
    return (0);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
