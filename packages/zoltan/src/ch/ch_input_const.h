// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __CH_INPUT_CONST_H
#define __CH_INPUT_CONST_H

#include <mpi.h>
#include "dr_const.h"
#include "dr_input_const.h"
#include "dr_compress_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

extern int chaco_input_graph(ZOLTAN_FILE*, char *, int **, int **, int *,
	   int *, float **, int *, float **);

extern int chaco_input_geom(ZOLTAN_FILE*, char *, int, int *, float **, float **,
			    float **);

extern int chaco_dist_graph(MPI_Comm, PARIO_INFO_PTR,
			    int, int *, int *, int **, int **,
			    int *, float **, int *, float **,
			    int *, float **, float **, float **, short **);
extern int chaco_input_assign(ZOLTAN_FILE*, char *, int, short *);


extern double read_val(ZOLTAN_FILE*, int *);
extern int read_int(ZOLTAN_FILE*, int *);

#ifndef TRUE
#define FALSE (0)
#define TRUE  (1)
#endif /* !TRUE */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
