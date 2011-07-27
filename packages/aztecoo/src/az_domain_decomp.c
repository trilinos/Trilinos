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

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


void AZ_domain_decomp(double x[], AZ_MATRIX *Amat, int options[],
                   int proc_config[], double params[],
		   struct context *context)


/*******************************************************************************

  Precondition 'x' using an overlapping domain decomposition method where a 
  solver specified by options[AZ_subdomain_solve] is used on the subdomains. 
  Note: if a factorization needs to be computed on the first iteration, this
  will be done and stored for future iterations.

  Author:          Lydie Prevost, SNL, 9222
  =======          Revised by R. Tuminaro (4/97), SNL, 9222

  Return code:     void
  ============

  Parameter list:
  ===============

  N_unpadded:      On input, number of rows in linear system (unpadded matrix) 
                   that will be factored (after adding values for overlapping).

  Nb_unpadded:     On input, number of block rows in linear system (unpadded) 
                   that will be factored (after adding values for overlapping).

  N_nz_unpadded:   On input, number of nonzeros in linear system (unpadded)
                   that will be factored (after adding values for overlapping).
             
  x:               On output, x[] is preconditioned by performing the subdomain
                   solve indicated by options[AZ_subdomain_solve].

  val    indx       
  bindx  rpntr:    On input, arrays containing matrix nonzeros (see manual). 
  cpntr  bpntr            

  options:         Determines specific solution method and other parameters.  In
                   this routine, we are concerned with options[AZ_overlap]:

                      == AZ_none: nonoverlapping domain decomposition
                      == AZ_diag: use rows corresponding to external variables 
                                  but only keep the diagonal for these rows.
                      == k      : Obtain rows that are a distance k away from
                                  rows owned by this processor.
                                  
  data_org:        Contains information on matrix data distribution and 
                   communication parameters (see manual).

*******************************************************************************/
{
  int N_unpadded, Nb_unpadded, N_nz_unpadded;
  double *x_pad = NULL, *x_reord = NULL, *ext_vals = NULL;
  int N_nz, N_nz_padded, nz_used;
  int mem_orig, mem_overlapped, mem_factor;
  int name, i, bandwidth;
  int *ordering = NULL;
  double condest;
/*
  double start_t;
*/
  int estimated_requirements;
  char str[80];
int *garbage;

  int N;
  int *padded_data_org = NULL, *bindx, *data_org;
  double *val;
  int *inv_ordering = NULL;
  int *map = NULL;
  AZ_MATRIX *A_overlapped = NULL;
  struct aztec_choices aztec_choices;


  /**************************** execution begins ******************************/
  data_org = Amat->data_org;
  bindx    = Amat->bindx;
  val      = Amat->val;
  N_unpadded = data_org[AZ_N_internal] + data_org[AZ_N_border];
  Nb_unpadded = data_org[AZ_N_int_blk]+data_org[AZ_N_bord_blk];
  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) 
     N_nz_unpadded = bindx[N_unpadded];
  else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
     N_nz_unpadded = (Amat->indx)[(Amat->bpntr)[Nb_unpadded]];
  else {
     if (Amat->N_nz < 0) 
        AZ_matfree_Nnzs(Amat);
     N_nz_unpadded = Amat->N_nz;
  }

  
  aztec_choices.options  = options;
  aztec_choices.params   = params;
  context->aztec_choices = &aztec_choices;
  context->proc_config   = proc_config;
  name                   = data_org[AZ_name];


  if ((options[AZ_pre_calc] >= AZ_reuse) && (context->Pmat_computed)) {
     N               = context->N;
     N_nz            = context->N_nz;
     A_overlapped    = context->A_overlapped;
     A_overlapped->data_org  = data_org;
     A_overlapped->matvec = Amat->matvec;
  }
  else {
     sprintf(str,"A_over %s",context->tag);
     A_overlapped = (AZ_MATRIX *) AZ_manage_memory(sizeof(AZ_MATRIX), 
                                                   AZ_ALLOC, name, str, &i);
     AZ_matrix_init(A_overlapped, 0);

     context->A_overlapped     = A_overlapped;
     A_overlapped->data_org    = data_org;
     A_overlapped->matvec      = Amat->matvec;
     A_overlapped->matrix_type = AZ_MSR_MATRIX;

     AZ_init_subdomain_solver(context);

     AZ_compute_matrix_size(Amat, options, N_nz_unpadded, N_unpadded, 
			 &N_nz_padded, data_org[AZ_N_external],
		 	 &(context->max_row), &N, &N_nz, params[AZ_ilut_fill], 
                         &(context->extra_fact_nz_per_row),
                         Nb_unpadded,&bandwidth);

     
        estimated_requirements = N_nz;
        if (options[AZ_overlap] != AZ_none && N_nz*2 > N_nz) {
          N_nz = N_nz*2;	/* check for overflow */
                          /* Add extra memory to N_nz. */
                          /* This extra memory is used */
                          /* as temporary space during */
                          /* overlapping calculation   */
        }

        /* Readjust N_nz by allocating auxilliary arrays and allocate */
        /* the MSR matrix to check that there is enough space.        */

        /* block off some space for map and padded_data_org in dd_overlap */

        garbage = (int *) AZ_allocate((AZ_send_list + 2*(N-N_unpadded) +10)*
                                      sizeof(int));
        AZ_hold_space(context, N);

   
        if (N_nz*((int) sizeof(double)) < N_nz) N_nz=N_nz/2; /* check for overflow */
        if (N_nz*((int) sizeof(double)) < N_nz) N_nz=N_nz/2; /* check for overflow */
        if (N_nz*((int) sizeof(double)) < N_nz) N_nz=N_nz/2; /* check for overflow */
        if (N_nz*((int) sizeof(double)) < N_nz) N_nz=N_nz/2; /* check for overflow */
        if (N_nz*((int) sizeof(double)) < N_nz) N_nz=N_nz/2; /* check for overflow */

        N_nz = AZ_adjust_N_nz_to_fit_memory(N_nz,
                                 context->N_large_int_arrays,
                                 context->N_large_dbl_arrays);
        context->N_nz_factors = N_nz;

        if (N_nz <= N_nz_unpadded) {
           AZ_printf_out("Error: Not enough space for domain decomposition\n");
           AZ_exit(1);
        }


        if (estimated_requirements > N_nz ) estimated_requirements = N_nz;

        /* allocate matrix via AZ_manage_memory() */

        sprintf(str,"bindx %s",context->tag);
        A_overlapped->bindx =(int    *) AZ_manage_memory(N_nz*sizeof(int),
                                                AZ_ALLOC, name, str, &i);
        sprintf(str,"val %s",context->tag);
        A_overlapped->val =(double *) AZ_manage_memory(N_nz*sizeof(double),
                                                AZ_ALLOC, name, str, &i);
        context->N_nz_allocated = N_nz;
        AZ_free_space_holder(context);
        AZ_free(garbage);

        /* convert to MSR if necessary */

        if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
          AZ_vb2msr(Nb_unpadded,val,Amat->indx,bindx,Amat->rpntr,Amat->cpntr,
		    Amat->bpntr, A_overlapped->val, A_overlapped->bindx);
        else if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
          for (i = 0 ; i < N_nz_unpadded; i++ ) {
             A_overlapped->bindx[i] = bindx[i];
             A_overlapped->val[i]   = val[i];
          }
        }
        else AZ_matfree_2_msr(Amat,A_overlapped->val,A_overlapped->bindx,N_nz);

        mem_orig = AZ_gsum_int(A_overlapped->bindx[N_unpadded],proc_config);

/*
        start_t = AZ_second();
*/
        AZ_pad_matrix(context, proc_config, N_unpadded, &N, 
                      &(context->map), &(context->padded_data_org), &N_nz, 
                      estimated_requirements);

/*
        if (proc_config[AZ_node] == 0)
           AZ_printf_out("matrix padding took %e seconds\n",AZ_second()-start_t);
*/


        mem_overlapped = AZ_gsum_int(A_overlapped->bindx[N],proc_config);
  
        if (options[AZ_reorder]) {
/*
           start_t = AZ_second();
*/
           AZ_find_MSR_ordering(A_overlapped->bindx,
                                &(context->ordering),N,
                                &(context->inv_ordering),name,context);
/*
           if (proc_config[AZ_node] == 0) 
              AZ_printf_out("took %e seconds to find ordering\n", AZ_second()-start_t);
*/
/*
           start_t = AZ_second();
*/
           AZ_mat_reorder(N,A_overlapped->bindx,A_overlapped->val,
                          context->ordering, context->inv_ordering);
/*
           if (proc_config[AZ_node] == 0) 
              AZ_printf_out("took %e seconds to reorder\n", AZ_second()-start_t);
*/
                /* NOTE: ordering is freed inside AZ_mat_reorder */
#ifdef AZ_COL_REORDER
           if (options[AZ_reorder]==2) {
             AZ_mat_colperm(N,A_overlapped->bindx,A_overlapped->val,
                        &(context->ordering), name, context );
           }
#endif

        }

        /* Do a factorization if needed.  */

/*
        start_t = AZ_second();
*/
        AZ_factor_subdomain(context, N, N_nz, &nz_used);

       if (options[AZ_output] > 0 && options[AZ_diagnostics]!=AZ_none) {
          AZ_printf_out("\n*********************************************************************\n");
	  condest = AZ_condest(N, context);
          AZ_printf_out("*****  Condition number estimate for subdomain preconditioner on PE %d = %.4e\n",
               proc_config[AZ_node], condest);
          AZ_printf_out("*********************************************************************\n");
        }


/*
        start_t        = AZ_second()-start_t;
        max_time = AZ_gmax_double(start_t,proc_config);
        min_time = AZ_gmin_double(start_t,proc_config);
        if (proc_config[AZ_node] == 0) 
           AZ_printf_out("time for subdomain solvers ranges from %e to %e\n",
                  min_time,max_time);
*/
  
        if ( A_overlapped->matrix_type == AZ_MSR_MATRIX)
           AZ_compress_msr(&(A_overlapped->bindx), &(A_overlapped->val),
                     context->N_nz_allocated, nz_used, name, context);


        context->N_nz = nz_used;
        context->N    = N;
        context->N_nz_allocated = nz_used;

        mem_factor = AZ_gsum_int(nz_used,proc_config);

        if (proc_config[AZ_node] == 0)
           AZ_print_header(options,mem_overlapped,mem_orig,mem_factor);

        if (options[AZ_overlap] >= 1) {
           sprintf(str,"x_pad %s",context->tag);
           context->x_pad  = (double *) AZ_manage_memory(N*sizeof(double),
                                                   AZ_ALLOC, name, str, &i);
           sprintf(str,"ext_vals %s",context->tag);
           context->ext_vals = (double *) AZ_manage_memory((N-N_unpadded)*
                                             sizeof(double), AZ_ALLOC, name, 
                                             str, &i);
        }
        if (options[AZ_reorder]) {
           sprintf(str,"x_reord %s",context->tag);
           context->x_reord = (double *) AZ_manage_memory(N*sizeof(double),
                                             AZ_ALLOC, name, str, &i);
        }

     }

  /* Solve L u = x where the solution u overwrites x */

    x_reord         = context->x_reord;
    inv_ordering    = context->inv_ordering;
    ordering        = context->ordering;
    x_pad           = context->x_pad;
    ext_vals        = context->ext_vals;
    padded_data_org = context->padded_data_org;
    map             = context->map;

   if (x_pad == NULL) x_pad = x;

   if (options[AZ_overlap] >= 1) {

      for (i = 0 ; i < N_unpadded ; i++) x_pad[i] = x[i];
      AZ_exchange_bdry(x_pad,padded_data_org, proc_config);
      for (i = 0 ; i < N-N_unpadded ; i++ ) 
         ext_vals[map[i]-N_unpadded] = x_pad[i+N_unpadded];
      for (i = 0 ; i < N-N_unpadded ; i++ ) x_pad[i + N_unpadded] = ext_vals[i];
   }
   else if (options[AZ_overlap] == AZ_diag) 
	AZ_exchange_bdry(x_pad,data_org, proc_config);

   if (x_reord == NULL) x_reord = x_pad;
   if (options[AZ_reorder]) {
      /* Apply row permutation to the right hand side */
      /* ((P'A P)Pi') Pi P'x = P'rhs,  b= P'rhs */
      for (i = 0 ; i < N ; i++ ) x_reord[inv_ordering[i]] = x_pad[i];
   }

   AZ_solve_subdomain(x_reord,N, context);

#ifdef AZ_COL_REORDER
   /* Apply column permutation to the solution   */
   if (options[AZ_reorder]==1){
      /* ((P'A P) P'sol = P'rhs   sol = P( P'sol)  */
      for (i = 0; i < N; i++) x_pad[i] = x_reord[inv_ordering[i]];
   }
   if (options[AZ_reorder]==2){
      /*
       * ((P'A P)Pi') Pi P'sol = P'rhs  sol = P Pi'( Pi P'sol)
       * Version 1:
       * for (i = 0; i < N; i++) pi_sol[i] = x_reord[ordering[i]];
       * for (j = 0; j < N; j++) x_pad[j] = pi_sol[inv_ordering[j]];
       * Version 2:
       */
      for (i = 0; i < N; i++) x_pad[i] = x_reord[ ordering[inv_ordering[i]] ];
   }

#else
   if (options[AZ_reorder])
      for (i = 0; i < N; i++) x_pad[i] = x_reord[inv_ordering[i]];
#endif

   AZ_combine_overlapped_values(options[AZ_type_overlap],padded_data_org, 
                             options, x_pad, map,ext_vals,name,proc_config);

   if (x_pad != x) 
     for (i = 0 ; i < N_unpadded ; i++ ) x[i] = x_pad[i];

} /* subdomain driver*/

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_print_header(int options[], int mem_overlapped,
                          int mem_orig, int mem_factor)
{
if ((options[AZ_overlap] < 1) && 
    (options[AZ_subdomain_solve] != AZ_ilut)) return;
 if ((options[AZ_output] != AZ_none ) && (options[AZ_output] != AZ_warnings) && (options[AZ_diagnostics]==AZ_all)){
      AZ_printf_out("\n\t\t*******************************************************\n");
      if (options[AZ_overlap] > 0) {
         AZ_printf_out("\t\t*****       Subdomain overlapping requires %.3e times\n", 
                ((double) mem_overlapped)/ ((double) mem_orig));
         AZ_printf_out("\t\t*****       the memory used for the nonoverlapped\n");
         AZ_printf_out("\t\t*****       subdomain matrix.\n");
      }
      if (options[AZ_subdomain_solve] == AZ_ilut) {
         AZ_printf_out("\t\t***** ilut: The ilut factors require %.3e times \n\t\t", 
                 ((double) mem_factor)/((double) mem_overlapped));
         AZ_printf_out("*****       the memory of the overlapped subdomain matrix.");
      }
      AZ_printf_out("\n\t\t*******************************************************\n");
      AZ_printf_out("\n");
   }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void AZ_find_MSR_ordering(int bindx2[],int **ordering,int N,
     int **inv_ordering, int name, struct context *context)

/*******************************************************************************
 
  Use a reverse cuthill McKee algorithm to find an ordering for the matrix.
 
  Author:          R. Tuminaro
 
  Return code:     void
  ============
 
  Parameter list:
  ===============
 
  bindx2:          On input, the nonzero sparsity pattern of the matrix
                   for which we will determine a new ordering. 
                   Note: bindx2 is changed in this routine, but then returned
                   to its original values before exiting.

  ordering:        On output, ordering[i] gives the new location of row i
                   in the reordered system.

  inv_ordering:    On output, inv_ordering[i] gives the location of row
*/
                   
{

   int i;
   int *mask;
   int root, nlvl, ccsize;
   int total = 0;
   char str[80];
  
   /* convert matrix to Fortran format */

   if (N==0) return;

   for (i = N+1 ; i < bindx2[N]; i++ ) bindx2[i]++;
   for (i = 0 ; i <= N ; i++ )         bindx2[i] -= N;

   /* initialize arrays for fnroot() and rcm() */

   sprintf(str,"inv_ordering %s",context->tag);
   *inv_ordering = (int *) AZ_manage_memory((N+1)*sizeof(int), AZ_ALLOC, name, 
                                            str,&i);
   sprintf(str,"ordering %s",context->tag);
   *ordering     = (int *) AZ_manage_memory((N+1)*sizeof(int), AZ_ALLOC, name,
                                            str,&i);
   mask          = (int *) AZ_allocate((N+1)*sizeof(int));
   if (mask == NULL) {
      AZ_printf_out("Not enough space for RCM reordering\n");
      AZ_exit(1);
   }

   for (i = 0 ; i < N ; i++ ) mask[i] = 1;
   root = 1;
   while (total != N ) {
      AZ_FNROOT_F77(&root,bindx2,&(bindx2[N+1]),mask, &nlvl, 
              &((*ordering)[total]), *inv_ordering);
      AZ_RCM_F77(&root,bindx2,&(bindx2[N+1]),mask,&((*ordering)[total]),
           &ccsize, *inv_ordering);

      if ( ccsize != N) {
         for (i = 0 ; i < ccsize ; i++) mask[(*ordering)[total+i]-1] = 0;
         for (i = 0 ; i < N ; i++ ) {
            if ( mask[i] == 1) break;
         }
         root = i+1;
      }
      total += ccsize;
      if (ccsize == 0) {
         AZ_printf_out("Error inside reordering\n");
         AZ_exit(1);
      }
   }

   /* convert matrix back to C format */

   for (i = 0 ; i <= N ; i++ ) bindx2[i] += N;
   for (i = N+1 ; i < bindx2[N]; i++ ) bindx2[i]--;

   /* convert ordering to C format */
   for (i = 0 ; i < N ; i++ ) (*ordering)[i]--;

   /* produce the inverse order */

   for (i = 0 ; i < N ; i++) (*inv_ordering)[(*ordering)[i]] = i;

   AZ_free(mask);
}

/* temp1 is some nz past the end of bindx[N] */
int AZ_pos_bin_search( int temp1, int bindx2[], int ordering[], int inv_ordering[],
                       double avg_nz_per_row, int N) 
{
#define COND(j) ((bindx2[(j)+1] > temp1))
#define COND_M(j) ((bindx2[(j)] <= temp1))
  /* define MY_ASSERT(cond,msg) do { if(!(cond)) { AZ_printf_out("ERROR: %s\n", msg); exit(1);   } } while(0) */
#define MY_ASSERT(cond,msg) ((void) (0))

  int irow = 0;
  int istart;
  int iend;
  int ifound = -1;
  int imid = -1;

  /* start binary search */
  {

    irow = (int) (((double)(temp1 - N))/avg_nz_per_row);

    istart = irow;
    iend = N-1;
    ifound = -1;
    imid = -1;

    /* forward search */
    if (COND(istart))
      {
        ifound = istart;
      }
    else
      {
        while (istart < iend - 1)
          {
            imid = (istart+iend) >> 1;
            if ( COND(imid) )
              {
                iend = imid;
              }
            else
              {
                istart = imid;
              }
          }

        if (COND(istart)) ifound = istart;
        if (COND(iend)) ifound = iend;
      }

    irow = ifound;

    /* backward search */
    istart = 0;
    iend = irow;
    ifound = -1;
    imid = -1;

    if (COND_M(iend))
      {
        ifound = iend;
      }
    else
      {
        while (istart < iend - 1)
          {
            imid = (istart+iend) >> 1;
            if (COND_M(imid))
              {
                istart = imid;
              }
            else 
              {
                iend = imid;
              }
          }

        if (COND_M(istart)) ifound = istart;
        if (COND_M(iend)) ifound = iend;
      }
    irow = ifound;
  }

  return( ordering[inv_ordering[irow]] + (temp1 - bindx2[irow]) );
}

/* index is some nz past the end of bindx[N] */
int AZ_pos( int index, int bindx[], int position[], int inv_ordering[],
            double avg_nz_per_row, int N) 
{
  int i = 0;

  i = (int) (((double)(index - N))/avg_nz_per_row);

  while ( bindx[i+1] <= index ) i++;
  while ( bindx[i] > index ) i--;

  return( position[inv_ordering[i]] + (index - bindx[i]) );
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_mat_reorder(int N, int bindx2_in[], double val2_in[], int ordering[],
	int inv_ordering[])
/*******************************************************************************
 
  Explicitly order a matrix.
 
  Author:          R. Tuminaro
 
  Return code:     void
  ============
 
  Parameter list:
  ===============
 
  bindx2:          On input, the nonzero sparsity pattern of the matrix
                   for which we will determine a new ordering. 
                   Note: bindx2 is changed in this routine, but then returned
                   to its original values before exiting.

  ordering:        On input, ordering[i] gives the new location of row i
                   in the reordered system.  Freed here.
                    
  inv_ordering:    On output, inv_ordering[i] gives the location of row
*/
{
  int irow = 0;
  int istart;
  int iend;
  int ifound = -1;
  int imid = -1;
  int endNZ = -1;

   int count, i, ordi;
   int mv_this_bindx, save_this_bindx;
   double mv_this_val, save_this_val;
   int current, temp1;
   double avg_nz_per_row;

   int * bindx2;  
   double * val2; 

   bindx2 = &bindx2_in[0];
   val2 = &val2_in[0];

   if (N == 0) return;
   avg_nz_per_row = ((double) (bindx2[N] - N))/((double) N);

   /* compute where row i's off-diagonals will be  */
   /* stored and put the result in position[i].    */

   count = N+1;
   for (i = 0 ; i < N ; i++ ) {
      ordi = ordering[i];
      ordering[i] = count;
      count += (bindx2[ordi+1] - bindx2[ordi]); 
   }

   /* move the offdiagonal elements in the rows */

   current = N+1;
   endNZ = bindx2[N];
   /*while (current < bindx2[N]) {*/
   while (current < endNZ) {
      mv_this_bindx   = bindx2[current];
      mv_this_val     = val2[current];
      temp1 = current;
      while ( mv_this_bindx >= 0 ) {
         /* temp1 = AZ_pos(temp1,bindx2,ordering,inv_ordering,avg_nz_per_row,N); */
         temp1 = AZ_pos_bin_search(temp1,bindx2,ordering,inv_ordering,avg_nz_per_row,N); 
         save_this_bindx = bindx2[temp1];
         save_this_val   = val2[temp1];
         bindx2[temp1] = -mv_this_bindx -1; /* encode this as a */
                                            /* negative value   */
         val2[temp1]   = mv_this_val;
         mv_this_bindx = save_this_bindx;
         mv_this_val   = save_this_val;
      }
      current++;
   }

   /* renumber the columns */

   /* for (i = N+1 ; i < bindx2[N] ; i++ ) bindx2[i]=inv_ordering[-bindx2[i]-1]; */
   for (i = N+1 ; i < endNZ ; i++ ) bindx2[i]=inv_ordering[-bindx2[i]-1];

   /* renumber the off-diagonal pointers */

   for (i = 0 ; i < N ; i++ ) bindx2[i] = ordering[i];

   for (i = 0 ; i < N; i++) ordering[i] = inv_ordering[i];
   AZ_sort(ordering, N, NULL, val2);

} /* end of AZ_mat_reorder */

void AZ_compute_matrix_size(AZ_MATRIX *Amat, int options[], int N_nz_unpadded,
	int N_unpadded, int *N_nz_padded, int N_external, int *max_row, int *N,
	int *N_nz, double fill, int *extra_fact_nz_per_row, int Nb_unpadded,
        int *bandwidth)
{

   int largest_padrow, extra_rows, extra_nonzeros, extra_factor_nonzeros;

   AZ_space_for_padded_matrix(options[AZ_overlap],N_nz_unpadded,N_unpadded,
                              &extra_rows, &extra_nonzeros,
                              N_external, &largest_padrow);

   *N           =  N_unpadded    + extra_rows;
   *N_nz        =  N_nz_unpadded + extra_nonzeros + 1;
   *N_nz_padded = *N_nz;
        
   *max_row = AZ_compute_max_nz_per_row(Amat,N_unpadded,Nb_unpadded,bandwidth);
   if (largest_padrow > *max_row ) *max_row = largest_padrow;

   AZ_space_for_factors(fill, *N_nz, *N, &extra_factor_nonzeros, options,
                        *bandwidth, *max_row);

   if (*N == 0) *extra_fact_nz_per_row = 0;
   else *extra_fact_nz_per_row = extra_factor_nonzeros/(2*(*N));
   *N_nz += extra_factor_nonzeros;
 
}

void AZ_pad_matrix(struct context *context, int proc_config[], 
   int N_unpadded, int *N, int **map, int **padded_data_org, 
   int *N_nz, int estimated_requirements)
{
   static int New_N_rows;
   int *data_org;
   int overlap;
   int i;
int *bindx;
double *val;
int count, start, end, j;

   data_org = context->A_overlapped->data_org;
   overlap  = context->aztec_choices->options[AZ_overlap];
   bindx    = context->A_overlapped->bindx;
   val      = context->A_overlapped->val;
   *map     = NULL; 
   *padded_data_org = data_org;

   if (overlap > 0) {
           *padded_data_org = data_org;
           New_N_rows = data_org[AZ_N_internal] + data_org[AZ_N_border];

           AZ_setup_dd_olap_msr(N_unpadded, &New_N_rows, bindx, val, overlap,
                           proc_config, padded_data_org,map, *N_nz, 
                           data_org[AZ_name], data_org, estimated_requirements,
			   context);

           if (New_N_rows > *N) {
              AZ_printf_out("Incorrectly estimated the overlap space reqirements.\n");
              AZ_printf_out("N_unpadded = %d, N_external = %d, overlap = %d\n",
		     N_unpadded, data_org[AZ_N_external], overlap);
              AZ_printf_out("Guess = %d, actual number of padded rows = %d\n",
                     *N, New_N_rows);
              AZ_printf_out("\n\nTry less overlapping and maybe we'll get it right\n");

              AZ_exit(1);
           }

           *N = New_N_rows;
    }
    else if (overlap == 0) {
       *N    = data_org[AZ_N_internal] + data_org[AZ_N_border];
       /* remove entries corresponding to external variables */

       count = bindx[0];
       start = count;
       for (i = 0 ; i < *N ; i++ ) {
          end = bindx[i+1];
          for (j = start ; j < end ; j++) {
             if ( bindx[j] < *N ) {
                bindx[count] = bindx[j];
                val[count++] = val[j];
             }
          }
          bindx[i+1] = count;
          start      = end;
       }

    }
    else { /* diagonal overlapping */

       *N = data_org[AZ_N_internal] + data_org[AZ_N_border];

       if (*N_nz < *N + data_org[AZ_N_external]) {
          AZ_printf_err("Not enough memory for diagonal overlapping\n");
          AZ_exit(1);
       }

       /* make room */

       count = data_org[AZ_N_external];
       for (i = bindx[*N]-1 ; i >= bindx[0] ; i-- ) {
          bindx[i+count] = bindx[i];
          val[i+count]   = val[i];
       }
       for (i = 0 ; i <= *N; i++) bindx[i] += count;
       for (i = (*N)+1 ; i <= *N + data_org[AZ_N_external]; i++) 
          bindx[i] = bindx[i-1];

       /* communicate diagonal */

       AZ_exchange_bdry(val, data_org, proc_config);

       *N = data_org[AZ_N_internal] + data_org[AZ_N_border] + 
                data_org[AZ_N_external];
    }
}

#ifdef abitmoregeneral
int AZ_pos( int index, int bindx[], int position[], int inv_ordering[],
         double avg_nz_per_row, int N, int first_row)
{
int i = 0;

   i = (int) (((double)(index - bindx[first_row]))/avg_nz_per_row);
   if (i >= N) i = N-1;
   i += first_row;
   while ( bindx[i+1] <= index ) i++;
   while ( bindx[i  ] >  index ) i--;
   return( position[inv_ordering[i-first_row]] + (index - bindx[i]) );
}

void AZ_mat2_reorder(int N, int first_row, int bindx2[], double val2[],
                     int ordering[], int inv_ordering[])
{
   int count, i, ordi;
   int mv_this_bindx, save_this_bindx;
   double mv_this_val, save_this_val;
   int current, temp1;
   double avg_nz_per_row;

   if (N == 0) return;
   avg_nz_per_row = ((double) (bindx2[first_row+N] - bindx2[first_row]))/
                                                ((double) N);

   /* compute where row i's off-diagonals will be  */
   /* stored and put the result in ordering[i].    */

   count = bindx2[first_row];
   for (i = 0 ; i < N ; i++ ) {
      ordi   = ordering[i] + first_row;
      ordering[i] = count;
      count += (  bindx2[ordi+1] - bindx2[ordi]);
   }

   /* move the offdiagonal elements in the rows */


   current = bindx2[first_row];
   while (current < bindx2[N + first_row]) {
      mv_this_bindx   = bindx2[current];
      mv_this_val     = val2[current];
      temp1 = current;
      while ( mv_this_bindx >= 0 ) {
         temp1 = AZ_pos( temp1,bindx2, ordering,inv_ordering,avg_nz_per_row,N,
                        first_row);
         save_this_bindx = bindx2[temp1];
         save_this_val   = val2[temp1];
         bindx2[temp1] = -mv_this_bindx -1; /* encode this as a */
                                            /* negative value   */
         val2[temp1]   = mv_this_val;
         mv_this_bindx = save_this_bindx;
         mv_this_val   = save_this_val;
      }
      current++;
   }
   /* renumber the columns */
   if (first_row == 0) {
      for (i = bindx2[0]; i < bindx2[N] ; i++ )
         bindx2[i]=inv_ordering[-bindx2[i]-1];
   }
   else {
      for (i = bindx2[first_row]; i < bindx2[first_row+N] ; i++ )
         bindx2[i]= -bindx2[i]-1;
   }

   /* renumber the off-diagonal pointers */

   for (i = 0 ; i < N ; i++ ) bindx2[i+first_row] = ordering[i];

   /* move around the diagonal elements */

   for (i = 0 ; i < N; i++) ordering[i] = inv_ordering[i];
   AZ_sort(ordering, N, NULL, &(val2[first_row]));

}
#endif
