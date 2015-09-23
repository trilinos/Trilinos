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

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id$"
;
#endif

/*******************************************************************************
 * MATRIX FREE  communication.
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
#include "example_specific.h"


void example_specific_comm_setup(int npoints, int *nx, int *ny, int proc, 
                             int nproc, struct exchange *comm_structure)

{
/*
 * Set up the data structure 'comm_structure'. In particular, set 
 *
 *   comm_structure->length_message[i] number of data points to exchange
 *                                     with the ith neighbor.
 *   comm_structure->neighbor_list[i]  proc id of ith neighbor (if the ith
 *                                     neighbor does not exist, = -1).
 *   comm_structure->buf               buffer space for grouping data
 *   comm_structure->send_list[i][j]   jth point index to send to ith neighbor
 *   comm_structure->recv_data[i]      data_buffer where information from the
 *                                     ith neighbor will be received.
 *
 * Note: i = LEFT, RIGHT, BOTTOM, TOP
 *
 * The data (npoints) is set so up as uniformly as possible on a processor 
 * array of size 'nproc X nproc'. 
 *
 * Parameters
 * ==========
 *
 * npoints          On input, number of data points to be distributed.
 *
 * nx,ny            On output, grid size of local grid.
 *
 * proc             On input, processor id of this processor.
 *
 * nproc            On input, nproc * nproc is the total number of processors.
 * 
 * comm_structure   On output, comm_structure is set (as described above) 
 *                  to do local data exchanges to update ghost variables.
 */



   int    px, py, N_augmented_procs,  j;
   int    send_coord;

  /* Convert the proc id's to processor array coordinates: (px, py) */
  /* Note: the processors are aligned on a nproc X nproc processor  */
  /*       array.                                                   */

  py = proc/nproc;
  px = proc%nproc;

  /* Compute the number of points per processor in each direction   */
  /* Note: some processors may contain 1 additional point in a      */
  /*       coordinate direction than others.                        */

  *nx = npoints/nproc;
  *ny = *nx;

  N_augmented_procs = npoints%nproc;

  if (px >= nproc-N_augmented_procs) *nx = *nx+1;
  if (py >= nproc-N_augmented_procs) *ny = *ny+1;


  /* allocate space to hold indices of data points to be */
  /* sent to neighbors and indices of data points to be  */
  /* received from neighbors.                            */

  comm_structure->send_list[LEFT  ] = (int *) malloc(*ny * sizeof(int));
  comm_structure->send_list[RIGHT ] = (int *) malloc(*ny * sizeof(int));
  comm_structure->send_list[BOTTOM] = (int *) malloc(*nx * sizeof(int));
  comm_structure->send_list[TOP   ] = (int *) malloc(*nx * sizeof(int));
  comm_structure->rcv_data[LEFT   ] = (double *) malloc(*ny * sizeof(double));
  comm_structure->rcv_data[RIGHT  ] = (double *) malloc(*ny * sizeof(double));
  comm_structure->rcv_data[BOTTOM ] = (double *) malloc(*nx * sizeof(double));
  comm_structure->rcv_data[TOP    ] = (double *) malloc(*nx * sizeof(double));

  if (comm_structure->rcv_data[TOP    ] == NULL) {
     printf("%d: Not enough memory to setup communication structure\n", proc);
     exit(1);
  }

  comm_structure->length_message[LEFT]   = 0;
  comm_structure->length_message[RIGHT]  = 0;
  comm_structure->length_message[TOP]    = 0;
  comm_structure->length_message[BOTTOM] = 0;
    
  /* get the neighbor numbers and indices */
  /* to exchange with each neighbor       */
    
  /* if there is a left neighbor*/

  comm_structure->neighbor_list[LEFT] = -1;
  if (px !=0) {
       comm_structure->neighbor_list[LEFT]  = proc-1;
       comm_structure->length_message[LEFT] = *ny;
       send_coord = 0;
       for ( j = 0; j< *ny; j++) {
             comm_structure->send_list[LEFT][j] = send_coord;
             send_coord += *nx;
       }
  }

  /* if right neighbor exists */

  comm_structure->neighbor_list[RIGHT] = -1;
  if (px != nproc -1) {
       comm_structure->neighbor_list[RIGHT]  = proc+1;
       comm_structure->length_message[RIGHT] = *ny;
       send_coord = *nx - 1;
       for ( j = 0; j< *ny; j++) {
          comm_structure->send_list[RIGHT][j] = send_coord;
          send_coord = send_coord + *nx;
       }
  }

  /* if bottom neighbor exists */

  comm_structure->neighbor_list[BOTTOM] = -1;
  if (py != 0) {
       comm_structure->neighbor_list[BOTTOM]  = proc-nproc;
       comm_structure->length_message[BOTTOM] = *nx;
       for ( j = 0; j< *nx; j++) 
          comm_structure->send_list[BOTTOM][j] = j;
  }

  /* if top neighbor exists */

  comm_structure->neighbor_list[TOP] = -1;
  if (py !=nproc -1) {
       comm_structure->neighbor_list[TOP] = proc+nproc;
       comm_structure->length_message[TOP] = *nx;
       send_coord = (*ny-1)*(*nx);
       for ( j = 0; j< *nx; j++) {
             comm_structure->send_list[TOP][j] = send_coord;
             send_coord++;
       }
  }
  comm_structure->buf   = (double *) malloc(2*((*nx)+(*ny)) * sizeof(double));

  if (comm_structure->buf  == NULL) {
     printf("%d: Not enough memory to setup communication structure\n",proc);
     exit(1);
  }
#ifdef AZTEC_MPI
  AZ_set_proc_config(comm_structure->proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(comm_structure->proc_config, AZ_NOT_MPI);
#endif


}


/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

void example_specific_exchange_data ( double *x, struct exchange *comm_structure)
{ 
/*
 * Exchange data between neighboring processors to update ghost
 * variables. Specifically, on output comm_structure->recv_data[i] holds
 * the appropriate data for each neighbor. This data is determined
 * by the structure 'comm_structure':
 *
 *     comm_structure->length_message[i]    number of data points to exchange
 *                                   with the ith neighbor.
 *     comm_structure->nieghbor_list[i]     proc id of ith neighbor (if the ith
 *                                   neighbor does not exist, = -1).
 *     comm_structure->buf                  buffer space for grouping data
 *     comm_structure->send_list[i][j]      index of jth point to send to ith neighbor
 *     comm_structure->recv_data[i]         data_buffer where information from the
 *                                   ith neighbor will be received.
 * Note: i = LEFT, RIGHT, BOTTOM, TOP
 */

int         *neighbor_list, *length_message; 
int         *send_list;
int         i, j, k, type, flag, count;
double      *buf; 
MPI_AZRequest request[4];
int         *proc_config;


    length_message = comm_structure->length_message;
    neighbor_list  = comm_structure->neighbor_list;
    buf            = (double *) comm_structure->buf;
    proc_config    = comm_structure->proc_config;
    type = 123;

    for (i = 0 ; i < 4 ; i++ ) {
       if (neighbor_list[i] != -1) 
          mdwrap_iread((void *) comm_structure->rcv_data[i], 
                        sizeof(double)*length_message[i], &(neighbor_list[i]), 
                        &type, &(request[i]));
    }   

    count = 0;
    for (i = 0 ; i < 4 ; i++ ) {
       if (neighbor_list[i] != -1) {
          send_list   =  comm_structure->send_list[i];
          k = count;
          for ( j =0;j<length_message[i];j++) buf[count++] = x[send_list[j]];
          mdwrap_write((void *) &(buf[k]), sizeof(double)*length_message[i],
		        neighbor_list[i], type, &flag);
       }
    }

    for ( i =0; i< 4; i++) {
       if (neighbor_list[i] != -1) {
          mdwrap_wait((void *) comm_structure->rcv_data[i], sizeof(double)*
                     length_message[i], &(neighbor_list[i]), &type, &flag,
                     &(request[i]));
       }
    }
}


/******************************************************************/
/******************************************************************/
/******************************************************************/
void example_specific_frees(struct exchange *comm_structure)
{
  /* free all fields of 'comm_structure' that we allocated inside */
  /* az_mat_free_com_set_up.                               */

  free(comm_structure->buf);
  free(comm_structure->rcv_data[TOP    ]);
  free(comm_structure->rcv_data[BOTTOM ]);
  free(comm_structure->rcv_data[RIGHT  ]);
  free(comm_structure->rcv_data[LEFT   ]);
  free(comm_structure->send_list[TOP   ]);
  free(comm_structure->send_list[BOTTOM]);
  free(comm_structure->send_list[RIGHT ]);
  free(comm_structure->send_list[LEFT  ]);
}
