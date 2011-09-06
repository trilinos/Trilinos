/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <float.h>
#include "az_aztec.h"

extern int AZ_sys_msg_type;
extern int AZ_using_fortran;


void AZ_read_update_funky(int *N_update, int *update[], int proc_config[],
                    int N, int chunk, int input_option)

/*******************************************************************************

  This routine initializes update[] to the global indices updated by this
  processor and initializes N_update to the total number of elements to be
  updated.

  If input_option == AZ_linear Do a linear partitioning of the chunks.
     Specifically, proc 0 is assigned the first floor( (N+P-1)/P ) chunks,
     processor 1 is assigned the next floor( (N+P-2)/P ) chunks, etc. where P =
     proc_config[AZ_N_procs].
  If input_option == AZ_file Processor 0 reads the file '.update'.  This file
     should contain nprocs lists.  Each list consists of a number telling how
     many global indices are in the list followed by a list of global indices.
     The first list is then sent to processor 'nprocs-1', the second list is
     sent to processor 'nprocs-2', etc.
  If input_option == AZ_box
     we do a box partitioning of the unknowns (see comments below).

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        On Output, number of unknowns updated by this processor.

  update:          On Output, list of unknowns updated by this processor in
                   ascending order.

  proc_config:     proc_config[AZ_node] is node number.
                   proc_config[AZ_N_procs] is the number of processors.

  N:               Total number of chunks to be distributed.

  chunk:           Size of each chunk to be treated as a single unit.
                   The unknowns contained in the kth chunk are given
                   by {k*chunk, k*chunk + 1, ..... , (k+1)*chunk - 1}
                   and 'N*chunk' is the total number of unknowns to be
                   distributed.

  input_option:    AZ_linear   ==> perform linear partitioning
                   AZ_file     ==> read partioning from file '.update'
                   AZ_box      ==> perform a box partitioning.

*******************************************************************************/

{

  /* local variables */

  int   t1, t2, i;
  int   ii, j;
  int   allocated, length;
  int   cflag;
  int   partner;
  int   proc_x, proc_y, proc_z;
  int   pts_x, pts_y, pts_z;
  int   total_pts_x, total_pts_y;
  int   px, py, pz, k;
  int   start_x, start_y, start_z;
  int   end_x, end_y, end_z;
  int   pt_number;
  int   count, check;
  int   proc, nprocs;
  int   type, type2;
  FILE *fp = NULL;
  MPI_AZRequest request;  /* Message handle */


  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  /*
   * Figure out which chunks should be assigned to this processor using a box
   * decomposition.  That is, it is assumed that all the chunks are ordered
   * naturally corresponding to an m x m x m box where m = N^(1/3).  Boxes of
   * chunks are assigned to processors.
   *
   * NOTE: it is assumed that nprocs = 2^power and that the number of chunks in
   * each direction is divisible by the number of processors in each direction.
   */

  if (input_option == AZ_box) {

    /* determine the number of processors in each direction */

    if (proc == 0) {
       (void) fprintf(stdout,"Input the dimensions of the processor cube\n\n");
       (void) fprintf(stdout,"Enter the number of processors along x axis>");
       (void) fflush(stdout);
       scanf("%d",&proc_x);
       (void) fprintf(stdout,"Enter the number of processors along y axis>");
       (void) fflush(stdout);
       scanf("%d",&proc_y);
       (void) fprintf(stdout,"Enter the number of processors along z axis>");
       (void) fflush(stdout);
       scanf("%d",&proc_z);

       (void) fprintf(stdout,"Input the grid dimensions\n\n");
       (void) fprintf(stdout,"Enter the number of grid points along x axis>");
       (void) fflush(stdout);
       scanf("%d",&pts_x);
       (void) fprintf(stdout,"Enter the number of grid points along y axis>");
       (void) fflush(stdout);
       scanf("%d",&pts_y);
       (void) fprintf(stdout,"Enter the number of grid points along z axis>");
       (void) fflush(stdout);
       scanf("%d",&pts_z);
    }
    AZ_broadcast((char *) &proc_x, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &proc_y, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &proc_z, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &pts_x , sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &pts_y , sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &pts_z , sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) NULL   , 0          , proc_config, AZ_SEND);
 
    total_pts_x = pts_x;
    total_pts_y = pts_y;


    if ( proc_x*proc_y*proc_z != nprocs) {
        if (proc == 0) {
          (void) fprintf(stdout,"Error: %d x %d x %d != %d ",
 			 proc_x, proc_y, proc_z, nprocs);
          (void) fprintf(stdout," (total number of processors)\n");
        }
	exit(1);
    }

    if ( pts_x * pts_y * pts_z != N ) {
        if (proc == 0) {
          (void) fprintf(stdout,"Error: %d x %d x %d != %d ",
 			 pts_x, pts_y, pts_z, N);
          (void) fprintf(stdout," (total number of grid points)\n");
        }
	exit(1);
    }
    if ( pts_x%proc_x != 0 ) {
        if (proc == 0) {
          (void) fprintf(stdout,"Error: grid points along x axis are not an ");
          (void) fprintf(stdout,"even multiple of processors\n");
	  (void) fprintf(stdout,"       along x axis.");
        }
	exit(1);
    }
    if ( pts_y%proc_y != 0 ) {
        if (proc == 0) {
          (void) fprintf(stdout,"Error: grid points along y axis is not an ");
          (void) fprintf(stdout,"even multiple of processors\n");
	  (void) fprintf(stdout,"       along y axis.");
        }
	exit(1);
    }
    if ( pts_z%proc_z != 0 ) {
        if (proc == 0) {
          (void) fprintf(stdout,"Error: grid points along z axis is not an ");
          (void) fprintf(stdout,"even multiple of processors\n");
	  (void) fprintf(stdout,"       along z axis.");
        }
	exit(1);
    }
    pts_x = pts_x/proc_x;
    pts_y = pts_y/proc_y;
    pts_z = pts_z/proc_z;

    /* compute the number of elements per processor in each direction */

    *N_update = pts_x * pts_y * pts_z * chunk;
    if (!AZ_using_fortran) 
       *update     = (int *) AZ_allocate((*N_update)*sizeof(int));

    /* compute the lower left corner and the upper right corner */

    px = proc % proc_x;
    pz = (proc-px) / proc_x;
    py = pz % proc_y;
    pz = (pz-py) / proc_y;

    start_x = px * pts_x;
    end_x   = start_x + pts_x;
    start_y = py * pts_y;
    end_y   = start_y + pts_y;
    start_z = pz * pts_z;
    end_z   = start_z + pts_z;

    /* set update[] */

    count = 0;
    for (k = start_z; k < end_z; k++ ) {
      for (j = start_y; j < end_y; j++ ) {
        for (i = start_x; i < end_x; i++ ) {
          for (ii = 0; ii < chunk; ii++ ) {
            pt_number = (i + j * total_pts_x + k * total_pts_x * total_pts_y) * 
                            chunk + ii;
            (*update)[count++] = pt_number;
          }
        }
      }
    }
  }

  else if (input_option == AZ_linear) {

    /*
     * Figure out which chunks should be assigned to this processor for linear
     * partitioning.  This means that processor 0 is assigned the chunks
     * approximately corresponding to 0, ... , N/nprocs and processor 1 is
     * approximately assigned the chunks 1+N/nprocs to 2*N/nprocs.
     */

    if (proc%2==1) {
     *N_update = 0; /* Odd processors get no update values */
     *update = NULL;
    }

    else {

    t1 = N/nprocs;
    t2 = N - t1 * nprocs;

    if ( proc >= t2) t2 += (proc * t1);
    else {
      t1++;
      t2    = proc*t1;
    }
    t1 *= 2; /* Double t1 to cover values from next processor */
    if ( proc+1 >= t2 && proc< t2) t1--; /* decrease t1 count by one if needed */
    *N_update = t1*chunk;
    t2   *= chunk;

    if (!AZ_using_fortran) 
       *update = (int *) AZ_allocate((*N_update+1)*sizeof(int));

    if (*update == NULL) {
      (void) fprintf (stderr, "Not enough space to allocate 'update'\n");
      exit(-1);
    }

    for (i = 0; i < *N_update; i++) (*update)[i] = i + t2;
    }
  }

  else if (input_option == AZ_file) {

    /* read the update elements from the file '.update' */

    type            = AZ_sys_msg_type;
    AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE)%AZ_NUM_MSGS +AZ_MSG_TYPE;
    type2           = AZ_sys_msg_type;
    AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE)%AZ_NUM_MSGS +AZ_MSG_TYPE;

    /*
     * Processor 0 reads file '.update' and distributes the lists to the other
     * processors.
     */

    t1 = 0;            /* total number of points distributed */
    if (proc == 0) {
      (void) printf("reading from file .update\n"); fflush(stdout);

      if ( (fp = fopen(".update", "r")) == NULL) {
        (void) fprintf(stderr, "ERROR: file '.update' not found\n");
        exit(-1);
      }

      if (!AZ_using_fortran) *update = 0;
      allocated   = 0;
      for (i = nprocs - 1; i >= 0; i-- ) {

        /* read in list length and send to processor i  */

        fscanf(fp, "%d", &length);
        t1 += length;
        if (i != 0)
          mdwrap_write((void *) &length, sizeof(int), i, type, &cflag);

        /*
         * If this is the last list, we allocate the right amount of space and
         * keep the list instead of sending it off
         */

        if (i == 0) {
          *N_update = length;
          allocated       = 0;
        }

        /* allocate enough space for list */

        if (length > allocated ) {
          if ((*update != NULL) && (!AZ_using_fortran)) AZ_free(*update);
          allocated = length + 1;

          if (!AZ_using_fortran)
            *update = (int *) AZ_allocate(allocated*sizeof(int));
          if (*update == NULL) {
            (void) fprintf(stderr,
                           "Not enough space to allocate 'update'\n");
            exit(-1);
          }
        }

        /* read a list and send it off to proc i (if not last list) */

        for (j = 0; j < length; j++ ) fscanf(fp, "%d", *update + j);
        if (i != 0)
          mdwrap_write((void *) *update, length*sizeof(int), i, type2, &cflag);
      }
      fclose(fp);

      if (t1 != N*chunk) {
        (void) fprintf(stderr,"AZ_read_update() found %d points in file\n", t1);
        (void) fprintf(stderr,"'.update' instead of the requested %d\n",
                       N*chunk);
        exit(-1);
      }
    }

    else {

      /* read the update list from processor 0 */

      partner = 0;
      mdwrap_iread((void *) N_update, sizeof(int), &partner, &type, &request);
      mdwrap_wait((void *) N_update, sizeof(int), &partner, &type, &cflag, &request);

      if (!AZ_using_fortran)
        *update = (int *) AZ_allocate((*N_update+1)*sizeof(int));
      if (*update == NULL)  {
        (void) fprintf(stderr, "Not enough space to allocate 'update'\n");
        exit(-1);
      }

      partner = 0;
      mdwrap_iread((void *) *update, *N_update * sizeof(int), &partner, &type2,
            &request);
      mdwrap_wait((void *) *update, *N_update * sizeof(int), &partner, &type2,
            &cflag, &request);
    }

    AZ_sort(*update, *N_update, NULL, NULL);

    /* check that element '0' is contained on 1 processor. That is,  */
    /* make sure the user has numbered from 0 to n-1 instead of from */
    /* 1 to n                                                        */
    check = 0;
    if ( (*N_update > 0) && ((*update)[0] == 0) ) check = 1;
    check = AZ_gsum_int(check, proc_config);
    if (check != 1) {
       if (proc == 0) {
          (void) fprintf(stderr,"Error: In AZ_read_update(), the '.update'");
          (void) fprintf(stderr,"file does not contain\n       one ");
          (void) fprintf(stderr,"occurance of row 0. Make sure that rows are");
          (void) fprintf(stderr," numbered\n       from 0 to n-1.\n");
       }
       exit(1);
    }


  }
  else {
    (void) fprintf(stderr,"Unknown input option (%d) in AZ_read_update()\n",
                   input_option);
    exit(1);
  }


} /* AZ_read_update */
