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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "az_aztec.h"
/* Begin Aztec 2.1 mheroux mod */
#ifdef IFPACK
#include "az_ifpack.h"
#include "az_blas_wrappers.h"
#endif
/* End Aztec 2.1 mheroux mod */

extern int az_iterate_id;

/*
 * 
 * These routines contain all the solver specific stuff for solving on subdomains.
 *
 * Note: See az_solve.c for an example of using these routines for an lu solver.
 *
 */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_factor_subdomain(struct context *context, int N, int N_nz,
	int *nz_used)
{
/****************************************************************************
  Given an overlapped subdomain matrix, factor it according to the
  chosen algorithm and store the result back in subdomain. Additionally, 
  store the number of nonzeros used in the factorization in nz_used.

  Notes:
    1) Matrix comes in as an MSR matrix.
    2) context contains several fields which need to be appropriately
       set. These fields are specific to the individual solvers.
    3) The factorization overwrites the matrix. However, different
       solvers will store the factorization in different formats.

  Author:          Ray Tuminaro, SNL, 9222 (3/98)

  Return code:     void
  ============

  Parameter list:
  ===============

  context        On input, context contains the matrix to be 
                   factored in context.A_overlapped (MSR format), 
                   On output, context contains the factored matrix
                   which is stored in a format specific to the solver and 
                   any additional parameters required by the backsolver.

  N                On input, the size of the linear system to be solved.

  N_nz             On input, the number of nonzero values in the matrix
                   to be factored.

  nz_used          On output, the number of nonzero values in the matrix
                   representing the factorization.

*******************************************************************************/

#ifdef HAVE_AZLU
	int ifail, N_nz_matrix, *rnr;
        double *fake_rhs, *aflag;
#endif
        int i, j, *bindx, *bpntr, *iw;
        double *cr, *unorm, *a, *val;
        int    *ind, *jnz, *ja, ifill;
        double dtemp = (context->aztec_choices->params)[AZ_omega];
        int    N_blk_rows, name = context->A_overlapped->data_org[AZ_name];
        char   str[80];

/* Begin Aztec 2.1 mheroux mod */
#ifdef IFPACK
	void *precon, *bmat;
	double rthresh, athresh;
	int N_int_blk, N_bord_blk, graph_fill;
#endif
/* End Aztec 2.1 mheroux mod */

        bindx = context->A_overlapped->bindx;

        *nz_used = bindx[N];

        switch(context->aztec_choices->options[AZ_subdomain_solve]) {
/* Begin Aztec 2.1 mheroux mod */
        case AZ_bilu_ifp:
#ifdef IFPACK
           if (N == 0) return;

           bindx = context->A_overlapped->bindx;
           val   = context->A_overlapped->val;

           /* for bilu(k) with k > 1 , figure out the new sparsity pattern */

           AZ_sort_msr(bindx, val, N);

           /* Let IFPACK handle fillin */
	      graph_fill = (context->aztec_choices->options)[AZ_graph_fill];
           (context->aztec_choices->options)[AZ_graph_fill] = 0;
           /* recover some space so that there will */
           /* be enough room to convert back to vbr */

           i = AZ_compress_msr(&(context->A_overlapped->bindx), 
                         &(context->A_overlapped->val), context->N_nz_allocated,
                         *nz_used, name, context);
           context->N_nz = *nz_used;
           context->N_nz_allocated = *nz_used;

           AZ_msr2vbr_mem_efficient(N, &(context->A_overlapped->bindx), 
                                 &(context->A_overlapped->val), 
                                 &(context->A_overlapped->cpntr), 
                                 &(context->A_overlapped->bpntr), 
                                 &(context->A_overlapped->indx), &N_blk_rows, 
                                 (context->A_overlapped->data_org)[AZ_name],
                                 context->tag,i);

	   context->A_overlapped->matrix_type = AZ_VBR_MATRIX;
   
	   /*ifp_initialize();*/
  
	   /* Create IFPACK encapsulation of Amat */

	   context->A_overlapped->rpntr = context->A_overlapped->cpntr;
	   N_int_blk = context->A_overlapped->data_org[AZ_N_int_blk];
	   N_bord_blk = context->A_overlapped->data_org[AZ_N_bord_blk];
	   context->A_overlapped->data_org[AZ_N_int_blk] = N_blk_rows;
	   context->A_overlapped->data_org[AZ_N_bord_blk] = 0;
	   (context->aztec_choices->options)[AZ_graph_fill] = graph_fill;

	   az2ifp_blockmatrix(&bmat, context->A_overlapped); 

	   context->A_overlapped->data_org[AZ_N_int_blk] = N_int_blk;
	   context->A_overlapped->data_org[AZ_N_bord_blk] = N_bord_blk;

	   rthresh =  (context->aztec_choices->params)[AZ_rthresh];
	   athresh =  (context->aztec_choices->params)[AZ_athresh];
           ifill = (context->aztec_choices->options)[AZ_graph_fill];
	   ifp_preconditioner(&precon, bmat, IFP_BILUK, (double) ifill, 0.0,
			    IFP_SVD, rthresh, athresh);
        if ((context->aztec_choices->options)[AZ_output]>0) {
              ifp_biluk_stats(precon);
        }
	   context->precon = precon;
           break;

/* End Aztec 2.1 mheroux mod */

#else
        AZ_perror("IFPACK not linked.  Must compile with -DIFPACK");
#endif
        case AZ_bilu:
           if (N == 0) return;

           bindx = context->A_overlapped->bindx;
           val   = context->A_overlapped->val;

           /* for bilu(k) with k > 1 , figure out the new sparsity pattern */

           AZ_sort_msr(bindx, val, N);
           ifill = (context->aztec_choices->options)[AZ_graph_fill];
           if (ifill > 0) {
              *nz_used = AZ_fill_sparsity_pattern(context, ifill, 
                                                  bindx, val, N);

           }
           /* recover some space so that there will */
           /* be enough room to convert back to vbr */

           i = AZ_compress_msr(&(context->A_overlapped->bindx), 
                         &(context->A_overlapped->val), context->N_nz_allocated,
                         *nz_used, name, context);
           context->N_nz = *nz_used;
           context->N_nz_allocated = *nz_used;

           AZ_msr2vbr_mem_efficient(N, &(context->A_overlapped->bindx), 
                                 &(context->A_overlapped->val), 
                                 &(context->A_overlapped->cpntr), 
                                 &(context->A_overlapped->bpntr), 
                                 &(context->A_overlapped->indx), &N_blk_rows, 
                                 (context->A_overlapped->data_org)[AZ_name],
                                 context->tag,i);

	   context->A_overlapped->matrix_type = AZ_VBR_MATRIX;
   
           bindx = context->A_overlapped->bindx;
           bpntr = context->A_overlapped->bpntr;
           val   = context->A_overlapped->val;

	   sprintf(str,"ipvt %s",context->tag);
           context->ipvt  = (int *) AZ_manage_memory((N+1)*sizeof(int),
                                    AZ_ALLOC, name, str, &i);
           sprintf(str,"dblock %s",context->tag);
           context->dblock= (int *) AZ_manage_memory((N_blk_rows+1)*
                                                 sizeof(int), AZ_ALLOC, name,
                                                 str, &i);

           context->N_blk_rows = N_blk_rows;

           /* set dblock to point to the diagonal block in each block row */

           for (i = 0 ; i < N_blk_rows ; i++ ) {
              for (j = bpntr[i] ; j < bpntr[i+1] ; j++ ) {
                 if (bindx[j] == i) context->dblock[i] = j;
              }
           }

           AZ_fact_bilu(N_blk_rows, context->A_overlapped, context->dblock,
                        context->ipvt);
           break;

	case AZ_ilut:
           cr = (double *) AZ_allocate((2*N+3+context->max_row)*sizeof(int)+
                                     (2*N+2+context->max_row)*sizeof(double));
           if (cr == NULL) AZ_perror("Out of space in ilut.\n");
           unorm = &(cr[N+2]);
           a     = &(unorm[N]);
           ind   = (int *) &(a[context->max_row]);
           jnz   = &(ind[N+3]);
           ja    = &(jnz[N]);
           sprintf(str,"iu %s",context->tag);
           context->iu    = (int *) AZ_manage_memory((N+1)*sizeof(int),
                                             AZ_ALLOC, name, str, &i);
           AZ_fact_ilut(&N, context->A_overlapped, a, ja,
                        (context->aztec_choices->params)[AZ_drop], 
                        context->extra_fact_nz_per_row, N_nz - bindx[N],
                        context->iu,cr,unorm,ind, nz_used, jnz,
                        (context->aztec_choices->params)[AZ_rthresh],
                        (context->aztec_choices->params)[AZ_athresh]);
           AZ_free(cr);
           break;
	case AZ_ilu:
           dtemp = 0.0;
	case AZ_rilu:
           if (N == 0) return;
           sprintf(str,"iu %s",context->tag);
           bindx = context->A_overlapped->bindx;
           val   = context->A_overlapped->val;

           /* for ilu(k) with k > 1 , figure out the new sparsity pattern */

           AZ_sort_msr(bindx, val, N);
           ifill = (context->aztec_choices->options)[AZ_graph_fill];
           if (ifill > 0) {
              *nz_used = AZ_fill_sparsity_pattern(context, ifill, 
                                                  bindx, val, N);
           }
           context->iu= (int *) AZ_manage_memory((N+1)*sizeof(int),AZ_ALLOC,
                                                    name, str, &i);
           iw = (int *) AZ_allocate((N+1)*sizeof(int));
           if (iw == NULL) AZ_perror("Out of space in ilu.\n");
           AZ_fact_rilu(N, nz_used, context->iu, iw, context->A_overlapped, 
                        dtemp,
                        (context->aztec_choices->params)[AZ_rthresh],
                        (context->aztec_choices->params)[AZ_athresh]);
           AZ_free(iw);
           break;
	case AZ_icc:
           sprintf(str,"iu %s",context->tag);
           bindx = context->A_overlapped->bindx;
           val   = context->A_overlapped->val;

           /* for ilu(k) with k > 1 , figure out the new sparsity pattern */

           AZ_sort_msr(bindx, val, N);
           ifill = (context->aztec_choices->options)[AZ_graph_fill];
           if (ifill > 0)
              *nz_used = AZ_fill_sparsity_pattern(context, ifill, 
                                                  bindx, val, N);

           AZ_fact_chol(context->A_overlapped->bindx,
                        context->A_overlapped->val,N,
                        (context->aztec_choices->params)[AZ_rthresh],
                        (context->aztec_choices->params)[AZ_athresh]);
           break;
	case AZ_lu:
#ifdef HAVE_AZLU
           if (N == 0) return;
           aflag = (double *) AZ_allocate(8*sizeof(double));
           rnr   = (int *) AZ_allocate(N_nz*sizeof(int));
           if (rnr == NULL) AZ_perror("Out of space in lu.\n");

           sprintf(str,"iflag %s",context->tag);
           context->iflag = (int *) AZ_manage_memory(10*sizeof(int), AZ_ALLOC,
                                                       name, str ,&i);
           sprintf(str,"ha %s",context->tag);
           context->ha = (int *) AZ_manage_memory(11*(N+1)*sizeof(int),
                                             AZ_ALLOC, name, str, &i);
           sprintf(str,"pivot %s",context->tag);
           context->pivot = (double *) AZ_manage_memory((N+1)*sizeof(double),
                                             AZ_ALLOC, name, str,&i);

           aflag[0] = 16.0;    aflag[2] = 1.0e8;   aflag[3] = 1.0e-12;   
           aflag[1] = (context->aztec_choices->params)[AZ_drop];

           /* set up flags for the sparse factorization solver */

           context->iflag[0] = 1;         context->iflag[1] = 2;
           context->iflag[2] = 1;         context->iflag[3] = 0;
           context->iflag[4] = 2;    
           /*    Note: if matrix is pos def, iflag[2] = 2 is cheaper */

           N_nz_matrix = bindx[N] - 1;

           AZ_msr2lu(N, context->A_overlapped, rnr);

           /* Mark bindx so we can see what was not used later */

           for (i =  N_nz_matrix ; i < N_nz ; i++) bindx[i] = -7;

           /* factor the matrix */ 

           if (N == 1) {
             context->A_overlapped->val[0]=1./context->A_overlapped->val[0];
           }
           else {
              context->N_nz_factors = N_nz;
              fake_rhs = (double *) AZ_allocate(N*sizeof(double));
              if (fake_rhs == NULL) {
                 AZ_printf_out("Not enough memory inside subdomain_solve\n");
              }
              for (i = 0 ; i < N ; i++ ) fake_rhs[i] = 0.0;
              AZ_fact_lu(fake_rhs, context->A_overlapped,aflag, 
                         context->pivot, rnr, context->ha, 
			 context->iflag, &N_nz_matrix,
                         &ifail, &(context->N_nz_factors),
                         &N, &N);

              (context->iflag)[4] = 3; 
              AZ_free(fake_rhs);

              /* find out what was not used by checking what was not touched */

              *nz_used = N_nz;
              for (i = N_nz_matrix; i < N_nz ; i++ ) {
                 if (bindx[i] != -7) *nz_used = i;
              }
              (*nz_used)++;
              context->N_nz_factors = *nz_used;
           }
           AZ_free(rnr);
           AZ_free(aflag);
#else
	   AZ_printf_err("AZ_lu unavailable: configure with --enable-aztecoo-azlu to make available\n");
	   exit(1);
#endif
           break;
        default:
           if (context->aztec_choices->options[AZ_subdomain_solve]
                  >= AZ_SOLVER_PARAMS) {
              AZ_printf_err("Unknown subdomain solver(%d)\n",
                   context->aztec_choices->options[AZ_subdomain_solve]);
              exit(1);
           }
        }      
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_free_space_holder(struct context *context)
{
/****************************************************************************
  This routine is used in conjunction with AZ_hold_space().
  Essentially, this routine deallocates memory allocated via
  AZ_hold_space(). The whole point of these two routines is to
  allocated all the space needed during the factorization process
  EXCLUDING all arrays whose size is related to the number of 
  nonzeros. Once this is done, we can determine how much space 
  there is left for the large arrays required for the factorization
  and split the remaining space amoung these large arrays. In this
  way LU routines where it is difficult to know the space requirements
  ahead of time can try to use as large an array as possible.

  Note: after factorization 'realloc' is used to reduce the array sizes.
  

  Author:          Ray Tuminaro, SNL, 9222 (3/98)

  Return code:     void
  ============

  Parameter list:
  ===============

  context          On input, context->aztec_choices->
		   options[AZ_subdomain_solve]
                   contains the preconditioner choice while
                   context->space_holder holds memory previously
                   allocated via AZ_hold_space().
                   On output, context->space_holder is deallocated.

                 
*******************************************************************************/
    int which = context->aztec_choices->options[AZ_subdomain_solve];

/* Begin Aztec 2.1 mheroux mod */
    if ( (which == AZ_ilut) || (which == AZ_lu ) || (which == AZ_bilu) ||
         (which == AZ_bilu_ifp) ||
         (which == AZ_rilu )|| (which == AZ_ilu) || (which == AZ_icc) )
       AZ_free(context->space_holder);
/* End Aztec 2.1 mheroux mod */
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_hold_space(struct context *context, int N) 
{
/****************************************************************************
  This routine is used in conjunction with AZ_free_space_holder().
  Essentially, this routine allocates memory while AZ_free_space_holder()
  deallocates.  The whole point of these two routines is to
  allocated all the space needed during the factorization process
  EXCLUDING all arrays whose size is related to the number of 
  nonzeros. Once this is done, we can determine how much space 
  there is left for the large arrays required for the factorization
  and split the remaining space amoung these large arrays. In this
  way LU routines where it is difficult to know the space requirements
  ahead of time can try to use as large an array as possible.

  Note: after factorization 'realloc' is used to reduce the array sizes.

  Author:          Ray Tuminaro, SNL, 9222 (3/98)

  Return code:     void
  ============

  Parameter list:
  ===============

  context          On input, context->aztec_choices->
	           options[AZ_subdomain_solve]
                   indicates the solver being used.
                   On output, context->space_holder points to a
                   block of memory which can hold all the 'nonlarge' arrays
                   required by this solver.
     
  N                On input, the size of the matrix to be allocated.

*******************************************************************************/


   switch(context->aztec_choices->options[AZ_subdomain_solve]) {
   case AZ_ilut:
      context->space_holder = (int *) AZ_allocate((2*N+2+context->max_row)*
                                                   sizeof(double) + sizeof(int)*
                                                   (3*N+8+context->max_row));
      if (context->space_holder  == NULL) AZ_perror("No space in ilut.\n");

      /*   Space for cr (N+2 doubles), unorm (N doubles), a (max_row doubles),*/
      /*   ind (N+3 ints), jnz (N ints), iu (N+1 ints), ja(max_row ints),     */
      /*   + 4 ints for manage memory header                                  */
      break;
#ifdef HAVE_AZLU
   case AZ_lu:
        context->space_holder = (int *) AZ_allocate((2*N+9)* sizeof(double) + 
                                           (11*(N+1)+22)*sizeof(int) );

      /*   Space for aflag (8 doubles), ifail (10 ints), ha (11*(N+1) ints), */
      /*   pivot (N+1 doubles), fake_rhs (N doubles)+12 ints for manage      */
      /*   memory header                                                     */
      /*   Note: Arrays of size N_nz (e.g. rnr) are not allocated by this    */
      /*   routine. Instead the subdomain field N_int_arrays is set.         */

        if (context->space_holder == NULL) AZ_perror("Out of space in lu.\n");
#else
	AZ_printf_err("AZ_lu unavailable: configure with --enable-aztecoo-azlu to make available\n");
	exit(1);
#endif
      break;
   case AZ_bilu:
/* Begin Aztec 2.1 mheroux mod */
   case AZ_bilu_ifp:
/* End Aztec 2.1 mheroux mod */
        context->space_holder= (int *) AZ_allocate((N+1)*sizeof(double));
        if (context->space_holder == NULL) AZ_perror("No space for bilu.\n");

        /* BILU is a little funny in that the maximum amount of memory */
        /* required does not occur during the factorization. Instead   */
        /* it occurs when converting MSR to VBR. At this point, an     */
        /* additional array of length N is needed.                     */
      break;

   case AZ_icc:
       context->space_holder= (int *) AZ_allocate((2*(N+1))*sizeof(int)+
                                                     (N+1)*sizeof(double));
       if (context->space_holder == NULL) AZ_perror("Out of space in ilu.\n");
       break;
   case AZ_ilu:
   case AZ_rilu:
       context->space_holder= (int *) AZ_allocate((2*(N+1)+4)*sizeof(int));

       /*   Space for iu (N+1 ints), iw (N+1 ints) + 4 ints for manage_memory */

       if (context->space_holder == NULL) AZ_perror("Out of space in ilu.\n");
       break;
    default:
      ;
   }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_init_subdomain_solver(struct context *context) 
{
/****************************************************************************
  Initialize the data structure 'context' to whatever values will
  be need by the particular routine.

  NOTE: The fields 'N_large_int_arrays' and 'N_large_dbl_arrays' must be 
  set to the number of arrays whose lengthes are equal to the number
  of nonzeros in the factorization. Most solvers will set these equal to
  one indicating that an msr bindx array is needed as well as an msr val
  array. However, the 'lu' solver for example requires an additional integer
  array to hold nonzeros and so N_large_int_arrays is set to 2. These values
  will be used when we try to allocate as much space as is needed to hold
  the factorization.

  Author:          Ray Tuminaro, SNL, 9222 (3/98)

  Return code:     void
  ============

  Parameter list:
  ===============

  context          On output, the various fields are set to solver specific
                   information that is needed in an initialization phase. 

  options          On input, is the user specified algorithm options given to 
                   AZ_solve() or AZ_iterate().

  params           On input, is the user specified algorithm options given to 
                   AZ_solve() or AZ_iterate().

*******************************************************************************/

   context->N_large_int_arrays = 1;
   context->N_large_dbl_arrays = 1;
   if (context->aztec_choices->options[AZ_subdomain_solve] == AZ_lu)
      context->N_large_int_arrays= 2;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_solve_subdomain(double x[],int N, struct context *context)
{
/****************************************************************************
  Given a vector 'x' representing the right hand side, solve the system
  using whatever subdomain solver is indicated by 'context->which'
  and whatever factorization information has already been computed.

  Author:          Ray Tuminaro, SNL, 9222 (3/98)

  Return code:     void
  ============

  Parameter list:
  ===============

  x                On input, the right hand side of the subdomain system that
                   is to be solved. 
                   On output, the solution of the subdomain system.

  N                On input, the size of the linear system to be solved.

  bindx2,val2      On input, matrix or factorization information to be used 
                   by the solver. For most schemes, this information is in
                   MSR format. However, the lu and bilu scheme would have
                   this information in another format.
                   Note: additional array information can be passed through
                   context.

  context          On input, the various fields are set to solver specific
                   information corresponding to algorithm parameters as
                   well as a previously done factorization.

*******************************************************************************/

double *val2;
int    *bindx2;
int N_blk_rows;
#ifdef HAVE_AZLU
int ifail;
#endif
int *sub_options, sub_proc_config[AZ_PROC_SIZE], *hold_data_org, *new_data_org;
double *sub_params, *sub_status;
AZ_MATRIX *sub_matrix;
AZ_PRECOND *sub_precond;
struct AZ_SCALING *sub_scaling;
#ifdef AZTEC_MPI
MPI_AZComm  *tptr;
#endif
double *y;
char label[80];
int  t1, t2, t3, i, t4, t5 = 0;

/* Begin Aztec 2.1 mheroux mod */
#ifdef IFPACK
  int ione = 1;
  void *precon;
#endif
/* End Aztec 2.1 mheroux mod */

   val2   = context->A_overlapped->val;
   bindx2 = context->A_overlapped->bindx;

   switch(context->aztec_choices->options[AZ_subdomain_solve]) {

/* Begin Aztec 2.1 mheroux mod */

   case AZ_bilu_ifp:
#ifdef IFPACK
     y = (double *) malloc (N * sizeof(double));
     DCOPY_F77(&N, x, &ione, y, &ione);
     precon = context->precon;
     ifp_apply(precon, N, 1, y, N, x, N);
     free((void *) y);
#endif
     break;

/* End Aztec 2.1 mheroux mod */

   case AZ_bilu:
      N_blk_rows = context->N_blk_rows;

      AZ_lower_triang_vbr_solve(N_blk_rows, context->A_overlapped->cpntr, 
                                context->A_overlapped->bpntr, 
				context->A_overlapped->indx,
                                bindx2, val2, x);

      AZ_upper_triang_vbr_solve(N_blk_rows, context->A_overlapped->cpntr,
                                context->A_overlapped->bpntr, 
				context->A_overlapped->indx, bindx2,
                                val2, x, context->ipvt, context->dblock);
      break;
   case AZ_ilut:
   case AZ_rilu:
   case AZ_ilu:
      AZ_lower_tsolve(x,N, val2, bindx2, context->iu, x ); 
      AZ_upper_tsolve( x, N, val2, bindx2, context->iu);
      break;
   case AZ_icc:
      AZ_lower_icc(bindx2,val2,N,x);
      AZ_upper_icc(bindx2,val2,N,x);
      break;
   case AZ_lu:
#ifdef HAVE_AZLU
      if (N == 0) return;
      else if (N== 1) {
         x[0] *= val2[0];
         ifail = 0;
      }
      else AZ_backsolve(val2, context->pivot,x, bindx2, 
	              context->ha, context->iflag, 
                      &ifail, &(context->N_nz_factors),
		      &N, &N);
#else
    AZ_printf_err("AZ_lu unavailable: configure with --enable-aztecoo-azlu to make available\n");
    exit(1);
#endif
      break;
   default: 
      if (context->aztec_choices->options[AZ_subdomain_solve]
                  >= AZ_SOLVER_PARAMS) {
         AZ_printf_out("ERROR: Unknown subdomain solver %d\n",
                context->aztec_choices->options[AZ_subdomain_solve]);
         exit(1);
       }
       else {
          /* better to put most of this in the factorization */

          AZ_recover_sol_params(context->aztec_choices->options[
			        AZ_subdomain_solve], &sub_options, 
				&sub_params, &sub_status, &sub_matrix, 
			        &sub_precond, &sub_scaling);
          t1 = sub_options[AZ_recursion_level];
          sub_options[AZ_recursion_level]++;

          t2 = sub_options[AZ_output];
          if (context->proc_config[AZ_node] != 0 ) 
             sub_options[AZ_output] = AZ_none;

          t3 = context->proc_config[AZ_MPI_Tag];

          /* fix data_org */

          hold_data_org = context->A_overlapped->data_org;
          new_data_org = (int *) AZ_allocate( sizeof(int) * AZ_send_list );
          if (new_data_org == NULL) {
             AZ_printf_out("Error: Not enough space for subdomain matrix\n");
             exit(1);
          }
          context->A_overlapped->data_org = new_data_org;
          context->A_overlapped->matvec = AZ_MSR_matvec_mult;
          new_data_org[AZ_matrix_type] = AZ_MSR_MATRIX;
          new_data_org[AZ_N_internal]  = N;
          new_data_org[AZ_N_border  ]  = 0;
          new_data_org[AZ_N_external]  = 0;
          new_data_org[AZ_N_int_blk ]  = N;
          new_data_org[AZ_N_bord_blk]  = 0;
          new_data_org[AZ_N_ext_blk ]  = 0;
          new_data_org[AZ_N_neigh   ]  = 0;
          new_data_org[AZ_total_send]  = 0;
          new_data_org[AZ_name      ]  = hold_data_org[AZ_name];
          new_data_org[AZ_internal_use]= 0;
          new_data_org[AZ_N_rows      ]= N;
          sub_precond->Pmat = context->A_overlapped;
          sub_precond->prec_function = AZ_precondition;
       
          sub_proc_config[AZ_node] = 0;
          sub_proc_config[AZ_N_procs] = 1;
#ifdef AZTEC_MPI
          tptr = AZ_get_comm(context->proc_config);
          AZ_set_comm(sub_proc_config, *tptr);
#endif

          sprintf(label,"y in ssolve%d", sub_options[AZ_recursion_level]);
          y = AZ_manage_memory((N+1)*sizeof(double),
                             AZ_ALLOC, AZ_SYS+az_iterate_id, label, &i);

          for (i = 0 ; i < N ; i++ ) y[i] = x[i];
          for (i = 0 ; i < N ; i++ ) x[i] = 0.0;

          t4 = sub_options[AZ_keep_info];
          sub_options[AZ_keep_info] = 1;

          if (context->aztec_choices->options[AZ_pre_calc] >= AZ_reuse) {
             t5 = sub_options[AZ_pre_calc];
             sub_options[AZ_pre_calc] = AZ_sys_reuse;
          }

          AZ_oldsolve(x, y,sub_options,sub_params, sub_status, sub_proc_config,
                       context->A_overlapped, sub_precond, sub_scaling);

          sub_options[AZ_keep_info] = t4;
          if (context->aztec_choices->options[AZ_pre_calc] == AZ_sys_reuse) 
             sub_options[AZ_pre_calc]  = t5;

          sub_options[AZ_recursion_level] = t1;
          sub_options[AZ_output] = t2;
          context->A_overlapped->data_org = hold_data_org;
          AZ_free(new_data_org);
          context->proc_config[AZ_MPI_Tag] = t3;
       }
   }
      
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_space_for_factors(double input_fill, int N_nz, int N, 
	int *extra_factor_nonzeros, int options[],int bandwidth,
        int  max_nz_per_row)
{
/****************************************************************************

  Compute the additional number of nonzeros required to do the factorization
  specified by options[AZ_subdomain_solve]. 
 
  Author:          Ray Tuminaro, SNL, 9222
 
  Return code:     void
  ============
 
  Parameter list:
  ===============
 
  input_fill:      On input, input_fill*N_nz is roughly the number of
                   nonzeros that will be allowed in the matrix factors
                   for ilut.

  N_nz:            On input, number of nonzeros in the padded matrix
                   that will be factored.
  
  N:               On input, the order of the padded matrix to be factored.

  extra_factor_nonzeros: 
                   On output, the additional space that will be added to
                   accommodate fill-in during the factorization.

  options:         On input, user specified input options.

 */

  int fill, i;
  double new_nz, nz_per_row, t;
  int    max_per_row, temp;
  
 
  if (options[AZ_subdomain_solve] == AZ_ilut) {
     input_fill -= 1.0;
     if (N == 0) *extra_factor_nonzeros = 0;
     else {

        new_nz = input_fill * ((double) N_nz);
        nz_per_row = new_nz/((double) N);
        fill = (int) floor( nz_per_row/2. + .5);
        t = ((double) N_nz)/ ((double) N);
        max_per_row = N - (int) ceil(t);
        if ( 2*fill > max_per_row) fill = max_per_row/2;

        *extra_factor_nonzeros =  2*N*fill + 1;
         while ( *extra_factor_nonzeros < 0) {
           fill--; *extra_factor_nonzeros =  2*N*fill + 1;
        }
     }
     temp = N*bandwidth;
     if (temp < 0) temp = *extra_factor_nonzeros;
     if (temp < *extra_factor_nonzeros) *extra_factor_nonzeros = temp;
  }
/* Begin Aztec 2.1 mheroux mod */
  else if ((options[AZ_subdomain_solve] == AZ_rilu) ||
           (options[AZ_subdomain_solve] == AZ_ilu ) ||
           (options[AZ_subdomain_solve] == AZ_icc ) ||
           (options[AZ_subdomain_solve] == AZ_bilu_ifp) ||
           (options[AZ_subdomain_solve] == AZ_bilu) ){
     fill = options[AZ_graph_fill];
/* End Aztec 2.1 mheroux mod */
     if (fill < 0) {
         AZ_printf_out("options[AZ_graph_fill] must be greater or equal to 0\n");
         AZ_printf_out("when using an incomplete factorization\n");
         exit(1);
     }
     if (fill == 0) *extra_factor_nonzeros = 3;
     else {
        temp = max_nz_per_row;
        for (i = 0 ; i < fill ; i++) { 
           temp *= max_nz_per_row;
           if (temp > bandwidth) break;
        }
        if (temp > bandwidth) temp = bandwidth;
        temp = temp*N;
        *extra_factor_nonzeros = temp - N_nz + 200;
     }
  }
  else if (options[AZ_subdomain_solve] == AZ_lu) {
     *extra_factor_nonzeros = N*bandwidth - N_nz + 200;
                             /* for small matrices y12m seems to need some    */
                            /* additional space. It might use a dense solver */
                            /* in some case???? who knows...                 */
  }
  else *extra_factor_nonzeros = 1;

  /* make sure things don't overflow */

  temp = 2*(N_nz + *extra_factor_nonzeros)*sizeof(double);
  if ( temp < 0 || *extra_factor_nonzeros < 0 ) {
      temp = 2;
      while ( temp < (2*temp) ) temp = 2*temp;
      *extra_factor_nonzeros = temp/(2*sizeof(double));
  }

}

void AZ_zero_out_context(struct context *context)
{

   context->iu     = NULL;
   context->iflag  = NULL;         
   context->ha     = NULL;
   context->ipvt   = NULL;
   context->dblock = NULL; 
   context->space_holder  = NULL;
   context->pivot  = NULL;
   context->A_overlapped  = NULL;
   context->aztec_choices = NULL;
   context->x_pad         = NULL; 
   context->ext_vals      = NULL; 
   context->x_reord       = NULL;
   context->padded_data_org = NULL;
   context->map           = NULL;
   context->inv_ordering  = NULL;
   context->tag           = NULL;
   context->proc_config   = NULL;

   context->extra_fact_nz_per_row = 0;
   context->N_large_int_arrays = 0;
   context->N_large_dbl_arrays = 0;
   context->N_nz_factors = 0;
   context->N_nz_matrix = 0; 
   context->N_blk_rows = 0;
   context->max_row = 0;
   context->N = 0;
   context->N_unpadded = 0;
   context->N_nz = 0; 
   context->Pmat_computed = 0;
}

