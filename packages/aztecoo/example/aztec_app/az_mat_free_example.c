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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
#include "example_specific.h"


void example_specific_matvec(double *x, double *y, AZ_MATRIX *Amat,
                        int proc_config[])
/******************************************************************************/
/*
 * Perform matrix vector operation:  
 * 
 *                    y = A x 
 *
 * where the matrix operator has the following stencil on
 * a 2D regular grid:
 *                                     -1
 *
 *                              -1      4     -1
 *    
 *                                     -1
 * Parameters:
 * =========
 * x                         On input, a vector which will be multiplied by
 *                           a matrix 'A'.
 *
 * y                         On output, the result of the matrix-vector 
 *                           product: A x.
 *
 * Amat                      On input, represents the matrix 'A'. Also
 *                           includes the values 'nx' and 'ny' which 
 *                           are the local grid size residing on the 
 *                           processor and the communication data structure
 *                           'comm_structure'. These are passed through using
 *                           AZ_get_matvec_data(Amat).
 * 
 */


{

   int nx, ny, px, py;
   int i, j, me, nproc;
   double temp;
   struct exchange *comm_structure;

   struct pass_data *pass_data;     /* Data passing structure. This user-     */
                                    /* defined data structure is used to pass */
                                    /* information through Aztec and back into*/
                                    /* the user's matrix-vector product.      */



  /*-------------------------------------------------------------------------*/

    pass_data      = (struct pass_data *) AZ_get_matvec_data(Amat);
    comm_structure = pass_data->comm_structure;
    nx = pass_data->nx;
    ny = pass_data->ny;

    nproc   = (double) sqrt( ((double) proc_config[AZ_N_procs]) + .01);
    py = (int) (proc_config[AZ_node]/nproc);
    px =  (int) (proc_config[AZ_node]%nproc);
    example_specific_exchange_data (x, comm_structure);

    for ( i = 0; i<nx; i++) {
      for ( j=0; j<ny;j++) {
         me = j*nx + i;

         temp = 4.*x[me];

         /* check in each direction to see if the current */
         /* grid point has a neighboring point in this    */
         /* direction (either on this processor or on a   */
         /* neighboring processor).                       */
         if      (i  != 0      ) temp -= x[me-1 ];
         else if (px != 0      ) temp -= comm_structure->rcv_data[LEFT][j];

         if      (i  != nx-1   ) temp -= x[me+1 ];
         else if (px != nproc-1) temp -= comm_structure->rcv_data[RIGHT][j];

         if      (j  != 0      ) temp -= x[me-nx];
         else if (py != 0      ) temp -= comm_structure->rcv_data[BOTTOM][i];

         if      (j  != ny-1   ) temp -= x[me+nx];
         else if (py != nproc-1) temp -= comm_structure->rcv_data[TOP][i];

         y[me] = temp;

      }
    }

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int  example_specific_diagonal_getrow( int columns[], double values[], 
	int row_lengths[], struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
        int requested_rows[], int allocated_space){
/*
 * Supply matrix diagonals for rows requested_rows[0 ... N_requested_rows-1].
 * Return this information in 'row_lengths, columns, values'.  If there is 
 * not enough space to complete this operation, return 0. Otherwise, return 1.
 *
 * Parameters     
 * ==========                                                        
 * Amat             On input, points to user's data containing matrix values.
 * N_requested_rows On input, number of rows for which nonzero are to be 
 *                  returned.                                  
 * requested_rows   On input, requested_rows[0...N_requested_rows-1] give the 
 *                  row indices of the rows for which nonzero values are 
 *                  returned.
 * row_lengths      On output, row_lengths[i] is the number of nonzeros in the 
 *                  row 'requested_rows[i]'  
 * columns,values   On output, columns[k] and values[k] contains the column
 *                  number and value of a matrix nonzero where all nonzeros for
 *                  requested_rows[i] appear before requested_rows[i+1]'s
 *                  nonzeros.  NOTE: Arrays are of size 'allocated_space'.
 * allocated_space  On input, indicates the space available in 'columns' and 
 *                  'values' for storing nonzeros. If more space is needed, 
 *                  return 0.
 */
   int i;

   if ( allocated_space < N_requested_rows ) return(0);

   for (i = 0; i < N_requested_rows; i++) {
      columns[i]     = requested_rows[i];
      values[i]      = 4.;
      row_lengths[i] = 1;
   }
   return(1);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int  simple_example_specific_getrow( int columns[], double values[], 
	int row_lengths[], struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
        int requested_rows[], int allocated_space){

/*
 * Supply local matrix (without ghost node columns) for rows given by
 * requested_rows[0 ... N_requested_rows-1].  Return this information in 
 * 'row_lengths, columns, values'.  If there is not enough space to complete 
 * this operation, return 0. Otherwise, return 1.
 *
 * Parameters     
 * ==========                                                        
 * Amat             On input, points to user's data containing matrix values.
 * N_requested_rows On input, number of rows for which nonzero are to be 
 *                  returned.                                  
 * requested_rows   On input, requested_rows[0...N_requested_rows-1] give the 
 *                  row indices of the rows for which nonzero values are 
 *                  returned.
 * row_lengths      On output, row_lengths[i] is the number of nonzeros in the 
 *                  row 'requested_rows[i]'  
 * columns,values   On output, columns[k] and values[k] contains the column
 *                  number and value of a matrix nonzero where all nonzeros for
 *                  requested_rows[i] appear before requested_rows[i+1]'s
 *                  nonzeros.  NOTE: Arrays are of size 'allocated_space'.
 * allocated_space  On input, indicates the space available in 'columns' and 
 *                  'values' for storing nonzeros. If more space is needed, 
 *                  return 0.
 */
   struct pass_data *pass_data;
   int    i,j,k, row, nz_ptr, old_ptr, nx, ny;

   pass_data = (struct pass_data *) AZ_get_getrow_data(Amat);
   nx        = pass_data->nx;
   ny        = pass_data->ny;

   nz_ptr = 0;
   for (k = 0; k < N_requested_rows; k++) {
      old_ptr = nz_ptr;
      row = requested_rows[k];
      i   = row%nx;
      j   = (row - i)/nx;

      if ( nz_ptr + 5 > allocated_space) return(0);

      values[nz_ptr   ] = 4.0;
      columns[nz_ptr++] = row;
      if      (i != 0    )    {values[nz_ptr]= -1.; columns[nz_ptr++]= row-1; }
      if      (i  != nx-1)    {values[nz_ptr]= -1.; columns[nz_ptr++]= row+1; }
      if      (j  != 0   )    {values[nz_ptr]= -1.; columns[nz_ptr++]= row-nx;}
      if      (j  != ny-1)    {values[nz_ptr]= -1.; columns[nz_ptr++]= row+nx;}

      row_lengths[k] = nz_ptr - old_ptr;
   }
   return(1);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int  example_specific_getrow( int columns[], double values[], 
	int row_lengths[], struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
        int requested_rows[], int allocated_space){

/*
 * Supply local matrix (with ghost node columns) for rows given by
 * requested_rows[0 ... N_requested_rows-1].  Return this information in 
 * 'row_lengths, columns, values'.  If there is not enough space to complete 
 * this operation, return 0. Otherwise, return 1.
 *
 * Parameters     
 * ==========                                                        
 * Amat             On input, points to user's data containing matrix values.
 * N_requested_rows On input, number of rows for which nonzero are to be 
 *                  returned.                                  
 * requested_rows   On input, requested_rows[0...N_requested_rows-1] give the 
 *                  row indices of the rows for which nonzero values are 
 *                  returned.
 * row_lengths      On output, row_lengths[i] is the number of nonzeros in the 
 *                  row 'requested_rows[i]'  
 * columns,values   On output, columns[k] and values[k] contains the column
 *                  number and value of a matrix nonzero where all nonzeros for
 *                  requested_rows[i] appear before requested_rows[i+1]'s
 *                  nonzeros.  NOTE: Arrays are of size 'allocated_space'.
 * allocated_space  On input, indicates the space available in 'columns' and 
 *                  'values' for storing nonzeros. If more space is needed, 
 *                  return 0.
 */
   struct pass_data *pass_data;
   int    i,j,k, row, nz_ptr, old_ptr, nx, ny, px, py, nproc;
   int    left, right, top, bottom;

   pass_data = (struct pass_data *) AZ_get_getrow_data(Amat);
   nx = pass_data->nx;
   ny = pass_data->ny;
   px = pass_data->px;
   py = pass_data->py;
   nproc = pass_data->nproc;

   left =  nx*ny;
   right = left;
   if (px != 0) right += ny;
   bottom = right;
   if (px != nproc-1) bottom += ny;
   top = bottom;
   if (py != 0) top += nx;

   nz_ptr = 0;
   for (k = 0; k < N_requested_rows; k++) {
      old_ptr = nz_ptr;
      row = requested_rows[k];
      i   = row%nx;
      j   = (row - i)/nx;

      if ( nz_ptr + 5 > allocated_space) return(0);

      values[nz_ptr]       = 4.0;
      columns[nz_ptr++] = row;
      if      (i != 0    )    {values[nz_ptr]= -1.; columns[nz_ptr++]= row-1; }
      else if (px != 0   )    {values[nz_ptr]= -1.; columns[nz_ptr++]= left+j;}

      if      (i  != nx-1)    {values[nz_ptr]= -1.; columns[nz_ptr++]= row+1; }
      else if (px != nproc-1) {values[nz_ptr]= -1.; columns[nz_ptr++]= right+j;}

      if      (j  != 0   )    {values[nz_ptr]= -1.; columns[nz_ptr++]= row-nx;}
      else if (py != 0   )    {values[nz_ptr]= -1.; columns[nz_ptr++]=bottom+i;}

      if      (j  != ny-1)    {values[nz_ptr]= -1.; columns[nz_ptr++]= row+nx;}
      else if (py != nproc-1) {values[nz_ptr]= -1.; columns[nz_ptr++]=top+i;  }

      row_lengths[k] = nz_ptr - old_ptr;
   }
   return(1);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int example_specific_comm_wrapper(double vec[], AZ_MATRIX *Amat)
/*
 * Update vec's ghost node via communication. Note: the length of vec is 
 * given by N_local + N_ghost where Amat was created via 
 *                 AZ_matrix_create(N_local);
 * and a 'getrow' function was supplied via
 *                 AZ_set_MATFREE_getrow(Amat, , , , N_ghost, );
 *
 * Parameters     
 * ==========                                                        
 * vec              On input, vec contains data. On output, ghost values
 *                  are updated.
 *
 * Amat             On input, points to user's data containing matrix values.
 *                  and communication information.
 */
{
   struct pass_data *pass_data;
   int    nx, ny, px, py, nproc, i, ptr;
   struct exchange *comm_structure;

   pass_data = (struct pass_data *) AZ_get_getrow_data(Amat);
   nx = pass_data->nx;
   ny = pass_data->ny;
   px = pass_data->px;
   py = pass_data->py;
   comm_structure = pass_data->comm_structure;
   nproc = pass_data->nproc;

   example_specific_exchange_data(vec, comm_structure);

   ptr = nx*ny;
   if (px != 0)
      for (i = 0; i < ny; i++) vec[ptr++] = comm_structure->rcv_data[LEFT][i];
   if (px != nproc-1)
      for (i = 0; i < ny; i++) vec[ptr++] = comm_structure->rcv_data[RIGHT][i];
   if (py != 0)
      for (i = 0; i < nx; i++) vec[ptr++] = comm_structure->rcv_data[BOTTOM][i];
   if (py != nproc-1)
      for (i = 0; i < nx; i++) vec[ptr++] = comm_structure->rcv_data[TOP][i];

   return(1);
}
