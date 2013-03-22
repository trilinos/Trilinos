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

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "dr_const.h"
#include "dr_loadbal_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************************************************************************/
void setup_fixed_obj(MESH_INFO_PTR mesh, int Num_Global_Parts)
{
int i, part;
int proc, nprocs;
FILE *fp;

  MPI_Comm_rank(MPI_COMM_WORLD, &proc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Initialize fixed elements in mesh */
  for (i = 0; i < mesh->num_elems; i++) mesh->elements[i].fixed_part = -1;

  switch (Test.Fixed_Objects) {
  case 1:
      /* Fix 100% of objects */
      for (i = 0; i < mesh->num_elems; i++) {
        mesh->elements[i].fixed_part = i % Num_Global_Parts;
      }
      break;
  case 2:
      /* Fix 50% of objects */
      for (i = 1; i < mesh->num_elems; i+=2) {
        mesh->elements[i].fixed_part = i % Num_Global_Parts;
      }
      break;
  case 3:
      /* Fix 10% of objects */
      for (i = 0; i < mesh->num_elems; i+=10) {
        mesh->elements[i].fixed_part = i % Num_Global_Parts;
      }
      break;
  case 4:
      /* Fix 100% of objects to partition 1 */
      for (i = 0; i < mesh->num_elems; i++) {
        mesh->elements[i].fixed_part = 1;
      }
      break;
  case 5:
      /* Fix 50% of objects to partition 1 */
      for (i = mesh->num_elems/2; i < mesh->num_elems; i++) {
        mesh->elements[i].fixed_part = 1;
      }
      break;
  case 6:
      /* Fix 10% of objects to partition 1 */
      for (i = 0; i < mesh->num_elems/10; i++) {
        mesh->elements[i].fixed_part = 1;
      }
      break;
  case 7:
      /* Fix one object on each proc. */
      for (i = 0; i < MIN(1, mesh->num_elems); i++) {
        mesh->elements[i].fixed_part = (proc%Num_Global_Parts);
      }
      break;
  case 8:
      /* Fix half the objects on half the procs, 25% overall */
      if (proc&1){
        for (i = 0; i < mesh->num_elems/2; i++) {
          mesh->elements[i].fixed_part = (proc%Num_Global_Parts);
        }
      }
      break;
  case 9:
      /* Fix first two objects on proc 0 only. */
      if (proc == 0){
        for (i = 0; i < MIN(2, mesh->num_elems); i++) {
          mesh->elements[i].fixed_part = i;
        }
      }
      break;
  case 10:
      /* Fix 0% of objects */
      break;

  case 99:
      /* Read from file. Assume serial execution for now. */
      if (proc == 0){
        /* TODO: Make filename a parameter. Hardcoded for now. */
        fp = fopen("fixed.dat", "r");
        if (fp == NULL)
           fprintf(stderr, "ERROR in opening file fixed.dat. No fixed vertices set.\n");
        else {
          while (1){
            if (feof(fp)) break; 
            /* read (i, part) for each fixed vertex */
            fscanf(fp, "%i%i\n", &i, &part);
            if ((part >= 0) && (part < Num_Global_Parts)){
              /* printf("Debug: setting fixed[%i] = %i\n", i, part); */
              mesh->elements[i].fixed_part = part;
            }
            else
              printf("Warning: Invalid part number %i ignored\n", part);
          }
        }
      }
      break;
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
