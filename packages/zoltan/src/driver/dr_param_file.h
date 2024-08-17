// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/*
   Code imported to Zoltan zdrive from

   zoltanParams.h

   prototypes for Zoltan parameters from a file, call Zoltan to set
   the parameters

   zoltanParams library

   Jim Teresco

   Department of Computer Science
   Williams College

   and

   Computer Science Research Institute
   Sandia National Laboratories

   Modification History:
   2004/02/02 JDT Created
   2004/03/10 JDT Removed obsolete hier_num_partitions function
*/

#ifndef __ZOLTANPARAMS_H
#define __ZOLTANPARAMS_H

#include <mpi.h>
#include "zoltan.h"

void zoltanParams_hier_free();
void zoltanParams_hier_set_num_levels(int levels);
void zoltanParams_hier_set_partition(int level, int partition);
void zoltanParams_hier_set_param(int level, char *param, char *value);
int zoltanParams_hier_get_num_levels();
int zoltanParams_hier_get_partition(int level);
void zoltanParams_hier_use_params(int level, struct Zoltan_Struct *zz, 
				  int *ierr);
void zoltanParams_set_comm(MPI_Comm thecomm);
void zoltanParams_hier_setup(struct Zoltan_Struct *zz);
void zoltanParams_read_file(struct Zoltan_Struct *lb, char *file, 
			    MPI_Comm thecomm);
#endif
