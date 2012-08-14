/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */
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
