/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

/* Data structure used for exchange information before performing */
/* a matrix-vector product in Aztec's matrix-free example.        */
/* NOTE:   NEIGHBOR(k) = (LEFT,RIGHT,BOTTOM,TOP)                  */

#define LEFT   0
#define RIGHT  1
#define BOTTOM 2
#define TOP    3

#define NO_PRECOND                  1
#define LS_PRECOND                  2
#define NS_PRECOND                  3
#define JAC_PRECOND                 4
#define NONOVERLAPDD_PRECOND        5
#define OVERLAPDD_PRECOND           6
#define USER_PRECOND                7


struct exchange {
   int neighbor_list[4];       /* [k]: proc id of NEIGHBOR(k)         */
   int length_message[4];      /* [k]: message length to exchange     */
                               /* with NEIGHBOR(k)                    */
   int    *send_list[4];       /* send_list[k][j]: index of jth point */
                               /*    to be sent to NEIGHBOR(k).       */
                            
   double *rcv_data[4];        /* rcv_data[i][j]: jth data point      */ 
                               /*    received from NEIGHBOR(k)        */
   double *buf;                /* buffer space for gathering messages */

   int    proc_config[AZ_PROC_SIZE];
};

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

struct pass_data {                  /* Data passing structure. This user-     */
   int nx, ny, px, py, nproc;       /* defined data structure is used to pass */
   struct exchange *comm_structure; /* information through Aztec and back into*/
};                                  /* the user's matrix-vector product.      */


/* -------------  external function declarations ------------------------- */

extern void example_specific_comm_setup (int N_coord, int *nx, int *ny, int ,
        int Nproc_coord, struct exchange *comm_structure);

extern void example_specific_matvec(double *x, double *y, AZ_MATRIX *Amat,
        int proc_config[]);

extern void example_specific_precond(double x[], int *, int *,
                     double *, AZ_MATRIX *, AZ_PRECOND *prec);

extern void example_specific_exchange_data(double *x, struct exchange
                                      *comm_structure);

extern int example_specific_comm_wrapper(double vec[], AZ_MATRIX *Amat);

extern void example_specific_frees(struct exchange *comm_structure);

extern int example_specific_read_input(int *options, int *proc_config,
                                        int *N_coord, int *Nproc_coord);

extern void example_specific_free(struct pass_data *pass_data);

extern int  example_specific_diagonal_getrow(int columns[], double values[],
        int row_lengths[], struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
        int requested_rows[], int allocated_space);

extern int  example_specific_getrow( int columns[], double values[],
        int row_lengths[], struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
        int requested_rows[], int allocated_space);

extern int example_specific_comm_wrapper(double vec[], AZ_MATRIX *Amat);

extern int  simple_example_specific_getrow( int columns[], double values[],
        int row_lengths[], struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
        int requested_rows[], int allocated_space);

