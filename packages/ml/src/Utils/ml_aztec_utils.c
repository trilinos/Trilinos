/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/************************************************************************/
/*          Utilities for Aztec/ML users                                */
/*----------------------------------------------------------------------*/
/* Author : Ray Tuminaro (SNL)                                          */
/************************************************************************/

#ifdef AZTEC
#ifdef ML_MPI
#define AZ_MPI
#endif
#include "ml_struct.h"
#include "ml_aztec_utils.h"

int warning_flag = 0;

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
int AZ_ML_Set_Amat(ML *ml_handle, int level, int isize, int osize, 
	AZ_MATRIX *Amat, int *proc_config)
{
/* Convert an Aztec matrix to an ML matrix and store the resulting ML
   matrix in the  level'th slot in ml_handle. 

   Note: This routine does not copy data. It simply associates wrapper
   functions to be used by ML to access data.

   Parameters
   ==========
   ml_handle             On input, ml_handle should have been created via
                         a prior call to ML_Create(). On output, the Aztec
                         matrix, Amat, is associated with the level'th 
                         matrix in ml_handle.

   level                 On input, level indicates where to store the 
                         information associated with Amat.

   isize, osize          On input, the length of the input and output vectors
                         (LOCALLY ON PROCESSOR NOT INCLUDING GHOST NODES) when
                         performing matrix-vector products.  Normally, these 
                         two lengths will be the same.

   Amat                  On input, an Aztec matrix that will be converted to 
                         an ML matrix.

   proc_config           On input, processor information (see Aztec guide).

*/
   struct aztec_context  *context;
   struct ML_CSR_MSRdata *msr_mat;
   struct ML_vbrdata     *vbr_mat;

   /* build Aztec context */

   context = (struct aztec_context *) ML_allocate(sizeof(struct aztec_context));
   context->Amat         = Amat;
   context->proc_config  = proc_config;
   context->comm         = ml_handle->comm;

   ML_Init_Amatrix(ml_handle, level,isize, osize, (void *) context);

   if (Amat->matrix_type == AZ_VBR_MATRIX) {
     vbr_mat = (struct ML_vbrdata *) AZ_allocate(sizeof(struct ML_vbrdata));
     vbr_mat->bindx       = Amat->bindx;
     vbr_mat->val         = Amat->val;
     vbr_mat->bpntr       = Amat->bpntr;
     vbr_mat->indx        = Amat->indx;
     vbr_mat->cpntr       = Amat->cpntr;
     vbr_mat->rpntr       = Amat->rpntr;
     context->getrowstuff = (void *) vbr_mat;
     ML_Set_Amatrix_Getrow(ml_handle,level,az_vbrgetrow_wrapper,az_comm_wrapper,
                           isize+(Amat->data_org)[AZ_N_external]);

     AZ_ML_set_vbrdiagonal(ml_handle,  level, Amat);
   }
   else if (Amat->matrix_type == AZ_MSR_MATRIX) {
     msr_mat = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct 
							      ML_CSR_MSRdata));
     msr_mat->columns     = Amat->bindx;
     msr_mat->values      = Amat->val;
     msr_mat->rowptr      = NULL;
     context->getrowstuff = (void *) msr_mat;
     ML_Set_Amatrix_Getrow(ml_handle,level,az_msrgetrow_wrapper,az_comm_wrapper,
                           isize+(Amat->data_org)[AZ_N_external]);
     ML_Set_Amatrix_Diag(  ml_handle, level, osize,   Amat->val);
   }
	 else if (Amat->matrix_type ==AZ_USER_MATRIX) {
		 context->getrowstuff = (void *)Amat->matvec_data;
     ML_Set_Amatrix_Getrow(ml_handle,level,az_usergetrow_wrapper,az_comm_wrapper,
													 isize+(Amat->data_org)[AZ_N_external]);
     AZ_ML_set_userdiagonal(ml_handle,  level, Amat);
	 }
   else {
      printf("Can only convert MSR, VBR or USER matrices\n");
      exit(1);
   }
   ML_Set_Amatrix_Matvec(ml_handle,  level, az_matvec_wrapper);
   ml_handle->Amat[level].data_destroy = AZ_ML_Clean;
   return(1);
}
int AZ_get_MSR_arrays(ML_Operator *Amat, int **bindx, double **val)
{
   struct aztec_context *context;
   struct ML_CSR_MSRdata *ptr = NULL;

   if (Amat->getrow->external == MSR_getrows) {
      ptr   = (struct ML_CSR_MSRdata *) Amat->data;
      *val   = ptr->values;
      *bindx = ptr->columns;
   }
   else if (Amat->getrow->external == az_msrgetrow_wrapper) {
      context = (struct aztec_context *) Amat->data;
      ptr = (struct ML_CSR_MSRdata *) context->getrowstuff;
      *val   = ptr->values;
      *bindx = ptr->columns;
   }
   else {
      *val   = NULL;
      *bindx = NULL;
/*
      printf("AZ_get_MSR_arrays: Not an msr matrix?\n");
      exit(1);
*/
   }
   return(1);
}

/***************************************************************************/
/*                     Wrapper for Aztec matvec                            */
/***************************************************************************/

int az_matvec_wrapper(void *data,  int in, double p[], int out, double ap[])
{
   struct aztec_context *temp;
   int      i, n,n2, *data_org;
   double   *p2;

   temp = (struct aztec_context *) data;
   data_org = temp->Amat->data_org;
   n        = data_org[AZ_N_internal] + data_org[AZ_N_border];
   n2       = n + data_org[AZ_N_external];
   p2       = (double *) AZ_allocate( (n2+1)*sizeof(double));
   for (i = 0; i < n; i++) p2[i] = p[i];
   temp->Amat->matvec(p2, ap, temp->Amat, temp->proc_config);
   for (i = 0; i < n; i++) p[i] = p2[i];
   AZ_free(p2);
   return(1);
}

/***************************************************************************/
/*                     Wrapper for Aztec communication                     */
/***************************************************************************/

int az_comm_wrapper(double vector[], void *data)
{
   struct aztec_context *temp;

   temp = (struct aztec_context *) data;

#ifndef AZTEC2_0
   AZ_exchange_bdry(vector, temp->Amat->data_org, temp->proc_config);
#else
   AZ_exchange_bdry(vector, temp->Amat->data_org);
#endif
   return 0;
}

/***************************************************************************/
/*                     Wrapper for Aztec MSR getrow                        */
/***************************************************************************/

int az_msrgetrow_wrapper(void *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   struct aztec_context *context;

   context = (struct aztec_context *) data;   

   return(MSR_getrows(context->getrowstuff, N_requested_rows, 
          requested_rows, allocated_space, columns, values, row_lengths) );
}

/***************************************************************************/
/*                     Wrapper for Aztec VBR getrow                        */
/***************************************************************************/

int az_vbrgetrow_wrapper(void *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   struct aztec_context *context;

   context = (struct aztec_context *) data;   

   return(VBR_cnst_blk_getrows(context->getrowstuff, N_requested_rows,
			       requested_rows,allocated_space,
			       columns, values, row_lengths) );
}

/***************************************************************************/
/*                     Wrapper for Aztec USER getrow                       */
/***************************************************************************/

int az_usergetrow_wrapper(void *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   struct aztec_context *context;
	 AZ_MATRIX *Amat;

   context = (struct aztec_context *) data;   

	 Amat=(AZ_MATRIX *)context->Amat;

   return(Amat->getrow(columns, values, row_lengths, Amat, N_requested_rows,
											 requested_rows,allocated_space));
}

/***************************************************************************/
/*            Memory deallocation for Aztec specific objects               */
/***************************************************************************/

void AZ_ML_Clean(void *data)
{
   struct aztec_context *context;

   context = (struct aztec_context *) data;
   if (context->Amat->matrix_type != AZ_USER_MATRIX) 
      ML_free(context->getrowstuff);
   ML_free(context);
} 


void AZ_ML_set_vbrdiagonal(ML *ml, int mesh_level, AZ_MATRIX *matrix)
{
/*
 *  Function to extract the diagonal entries from a VBR matrix and pass them
 *  to the ML object.  Author == Ray Tuminaro.
 */

  int i, j, k, m, off, start, end, fixed_leng, num_blks, blk_size;
  double *diagonal;

  fixed_leng = matrix->data_org[AZ_N_internal] + matrix->data_org[AZ_N_border];
  diagonal = (double *) ML_allocate( (fixed_leng+1)*sizeof(double));
  num_blks = matrix->data_org[AZ_N_int_blk] + matrix->data_org[AZ_N_bord_blk];
  for (i=0, k=0; k < num_blks; k++) {
     start = matrix->bpntr[k];
     end   = matrix->bpntr[k+1] - 1;
     for (j=start; j <= end; j++) {
        if ( matrix->bindx[j] == k ) break;
     }
     blk_size =  matrix->rpntr[k+1]-matrix->rpntr[k];
     for (m=0, off=0; m < blk_size; m++) {
        diagonal[i++] = matrix->val[matrix->indx[j] + off];
        off += blk_size + 1;
     }
  }
  ML_Set_Amatrix_Diag( ml, mesh_level, fixed_leng, diagonal );
  ML_free(diagonal);
}

void AZ_ML_set_userdiagonal(ML *ml, int mesh_level, AZ_MATRIX *matrix)
{
/*
 *  Function to extract the diagonal entries from a USER matrix and pass them
 *  to the ML object.  Author == Dawn Chamberlain.
 */

  int i, tmp, loc, fixed_leng, row_len, *cols, max_nnz_per_row=500;
  double *diagonal, *vals;

  fixed_leng = matrix->data_org[AZ_N_internal] + matrix->data_org[AZ_N_border];
  diagonal = (double *) ML_allocate( (fixed_leng)*sizeof(double));

	cols = (int *) malloc(max_nnz_per_row*sizeof(int));
	vals = (double *)malloc(max_nnz_per_row*sizeof(double));
	if (vals == NULL) {
		printf("AZ_ML_set_userdiagonal: memory allocation error\n");
		exit(-1);
	}


	for (i=0; i < fixed_leng; i++) {
		tmp = matrix->getrow(cols, vals, &row_len, matrix, 1, 
												 &i, max_nnz_per_row);
		while (tmp == 0) {
			free(cols);
			free(vals);
			max_nnz_per_row=max_nnz_per_row*2+1;
			cols=(int *)malloc(max_nnz_per_row*sizeof(int));
			vals=(double *)malloc(max_nnz_per_row*sizeof(double));
			tmp = matrix->getrow(cols, vals, &row_len, matrix, 1, 
													 &i, max_nnz_per_row);
		}
		loc=0;
		while ((loc<row_len) && (cols[loc] != i))
			loc++;

		if (loc == row_len)
			diagonal[i]=0.0;
		else
			diagonal[i]=vals[loc];
	}

		
  ML_Set_Amatrix_Diag( ml, mesh_level, fixed_leng, diagonal );
  ML_free(diagonal);
	free(cols);
	free(vals);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void ML_precondition(double ff[], int options[], int proc_config[],
                     double params[], AZ_MATRIX *mat, AZ_PRECOND *prec)
{
/*
 * Preconditioning wrapper function to be called by Aztec when using 
 * ML as a preconditioner.
 */

  int         i = 0, lenf;
  double      *ffout;
  static      int message = 0;
#ifdef ML_TIMING
  double      t0;
#endif

  ML          *ml;

#ifdef AZ_ver2_1_0_3
  ml    = (ML *) AZ_get_precond_data(prec);
#else
  ml    = (ML *) prec->ml_ptr;
#endif
#ifdef ML_TIMING
  t0 = GetClock();
#endif
  if (message == 0) {
     message  = 1;
     if (   (options[AZ_solver] != AZ_fixed_pt) &&
            (options[AZ_solver] != AZ_GMRESR) &&
            (warning_flag == 1) &&
            (ml->comm->ML_mypid == 0) ) {
        printf("Warning:Using a Krylov method to precondition a "); 
        printf("Krylov\n");
        printf("\tmethod has 'some' risk (as the preconditioner\n"); 
        printf("\tmight change from iteration to iteration).\n");
        printf("\tSetting options[AZ_solver] = AZ_GMRESR invokes an\n");
        printf("\tunsupported solver intended to handle changing \n");
        printf("\tpreconditioners or ML_Iterate() can be used to run\n");
        printf("\tthe multilevel method.\n\n"); 
     }
  }
  lenf  = ml->SingleLevel[ml->ML_finest_level].Amat->outvec_leng;

  /* then apply a two level preconditioning */

  ffout = (double*) malloc(lenf * sizeof(double));
  if (ml->ML_scheme == ML_MGFULLV) ML_Solve_MGFull( ml, ff, ffout );
  else ML_Solve_MGV( ml, ff, ffout );
  for (i = 0; i < lenf; i++) ff[i] = ffout[i];
  free(ffout);
#ifdef ML_TIMING
  ml->timing->precond_apply_time += (GetClock() - t0);
#endif
}

/*****************************************************************************/
/*****************************************************************************/

void AZ_set_ML_preconditioner(AZ_PRECOND **Precond, AZ_MATRIX *Amat, 
                              ML *ml_handle, int options[])
{
   char str[80], coarsest[160], finest[160];
   int  i;

   if (*Precond != NULL) {
      printf("AZ_set_ML_preconditioner: *Precond is not NULL. Is there already a preconditioner?\n");
      printf("\t\tIf so, use AZ_precond_destroy to remove. Otherwise, set to NULL before\n");
      printf("\t\tinvoking AZ_set_ML_preconditioner().\n");
      exit(1);
   }
#ifdef AZ_ver2_1_0_3
   *Precond = AZ_precond_create(Amat, ML_precondition, ml_handle);
   /*  take finest and coarsest grid smoothers and stick in string */
   /*  along with the number of levels                             */
   i = ml_handle->ML_finest_level;
   finest[0]   = '\0';
   coarsest[0] = '\0';
   if (i != -1) {
      if (ml_handle->pre_smoother[i].ML_id != ML_EMPTY) 
         sprintf(finest, "%s", ml_handle->pre_smoother[i].label);
      if (ml_handle->post_smoother[i].ML_id != ML_EMPTY) 
         sprintf(finest, "%s/%s", finest,ml_handle->post_smoother[i].label);

      if (i != ml_handle->ML_coarsest_level) {
         i = ml_handle->ML_coarsest_level;
         if (ml_handle->csolve[i].ML_id != ML_EMPTY) 
            sprintf(coarsest, "%s", ml_handle->csolve[i].label);
         else {
            if (ml_handle->pre_smoother[i].ML_id != ML_EMPTY) 
            sprintf(coarsest, "%s", ml_handle->pre_smoother[i].label);
            if (ml_handle->post_smoother[i].ML_id != ML_EMPTY) 
            sprintf(coarsest, "%s/%s", finest,ml_handle->post_smoother[i].label);
         }
      }
      sprintf(str,"%d level MG ( %s, %s)", ml_handle->ML_num_actual_levels, finest, coarsest);
   }
   else sprintf(str,"multilevel");
   AZ_set_precond_print_string(*Precond, str);
#else
   *Precond = AZ_precond_create(Amat, ML_precondition);
#endif
   options[AZ_precond]    = AZ_user_precond;
   (*Precond)->ml_ptr = (void *) ml_handle;
}

/*****************************************************************************/
/*****************************************************************************/
int AZ_block_MSR(int **param_bindx, double **param_val,
                 int N_update, int num_PDE_eqns, int *update)
{
/*
 * Take an MSR matrix and stick zeros into it so that it corresponds
 * to a block matrix where each block is num_PDE_eqns x num_PDE_eqns.
 * This is useful inside ML when we do Amalgmation for coarsening.
 */

   int allocated, nblocks, *block_list, i, j, k;
   int offset, current, block;
   int *bindx, *newbindx;
   double *val, *newval;

   bindx = *param_bindx;
   val   = *param_val;

   allocated  = (int     ) (((double)(bindx[N_update]+5))*1.2);
   block_list = (int    *) AZ_allocate( 300*sizeof(int));

   newbindx   = (int    *) AZ_allocate( allocated*sizeof(int));
   newval     = (double *) AZ_allocate( allocated*sizeof(double));
   *param_bindx = newbindx;
   *param_val   = newval;

   if (newval == NULL) {printf("out of space\n"); exit(1); }

   /* Get the diagonal */

   for (i = 0; i < N_update; i++) newval[i] = val[i];
   for (i = 0; i < N_update; i++) newbindx[i] = bindx[i+1] - bindx[i];

   offset = bindx[0];
   current = bindx[0];
   newbindx[0] = bindx[0];
   AZ_sort_msr(bindx,val,N_update);

   for (i = 0; i < N_update; i++) {

      /* make a list of all the block columns in this row */

      nblocks = 0;
      block_list[nblocks++] = update[i]/num_PDE_eqns;
      for (j = bindx[i]; j < bindx[i+1]; j++) {
         block = bindx[j]/num_PDE_eqns;
         if ((block != block_list[0]) && (block != block_list[nblocks-1])) {
            block_list[nblocks++] = block;
         }
      }

      /* sort the block and make sure that within    */
      /* each block we have each point column entry. */

      AZ_sort(block_list, nblocks, NULL, NULL);
      for (j = 0; j < nblocks; j++) {
         for (k = 0; k < num_PDE_eqns; k++) {
            if (block_list[j]*num_PDE_eqns+k != update[i]) {
               if (offset >= allocated)
                  pr_error("ML_block_MSR: Did not allocate enough space\n");

               newbindx[offset] = block_list[j]*num_PDE_eqns+k;
               if ( (current < bindx[i+1]) &&
                    (bindx[current] == newbindx[offset]) ) {
                  newval[offset++] = val[current++];
               }
               else { newval[offset++] = 0.;  }
            }
         }
      }
      newbindx[i+1] = offset;
   }
   return 0;
}


#ifndef AZTEC2_0

/*****************************************************************************/
/*****************************************************************************/

void ML_Gen_SmootherAztec(ML *ml_handle, int level, int options[], 
	double params[], int proc_config[], double status[],
	int N_iterations, int pre_or_post,
        void (*prec_function)(double *, int *, int *, double *,
                              struct AZ_MATRIX_STRUCT  *,
                              struct AZ_PREC_STRUCT *))
{
   struct aztec_context  *context, *orig_context;
   struct ML_CSR_MSRdata *msr_mat;
   ML_Operator           *ML_Amat;
   AZ_MATRIX             *AZ_Amat = NULL, *orig_Amat;
   AZ_PRECOND            *AZ_Prec = NULL;
   int                   *options_copy = NULL, i, *data_org = NULL;
   double                *params_copy = NULL;
   int                   save_old_values[6], *orig_data_org, global_mat_flag=0;
   static int            matrix_name = 7911;
   void           *data;
   ML_Operator    *op;
   int            osize, *row_ptr, space, getrow_flag, flag, *cols, nz_ptr;
   int            length, zero_flag, j, offset, nrows, *sub_proc_config;
   double         *vals, dsize, di;
   ML_Matrix_DCSR *csr_mat, *csr2_mat;
   #ifdef ML_MPI
      MPI_AZComm *tptr;
   #endif
#ifdef ML_EXPERIMENT
   /* invoke Aztec once to build preconditioner. */
   int size;
   double *xxx, *rhss;
#endif

   /* build Aztec context */

   if (options[AZ_scaling] != AZ_none) {
      if (ml_handle->comm->ML_mypid == 0) {
         printf("ML_Gen_SmootherAztec: Can not use Aztec scaling\n");
      }
      exit(1);
   }
   options_copy = (int    *) AZ_allocate(sizeof(int   )*AZ_OPTIONS_SIZE);
   params_copy  = (double *) AZ_allocate(sizeof(double)*AZ_PARAMS_SIZE);
   for (i = 0 ; i < AZ_OPTIONS_SIZE; i++) options_copy[i] = options[i];
   for (i = 0 ; i < AZ_PARAMS_SIZE ; i++) params_copy[i]  = params[i];
   options_copy[AZ_output]    = AZ_none;
   options_copy[AZ_recursion_level] = options[AZ_recursion_level] + 1;
   options_copy[AZ_keep_info] = 1;
   if (N_iterations != AZ_ONLY_PRECONDITIONER) warning_flag = 1;
   if (options_copy[AZ_max_iter] == -42) global_mat_flag = 1;
   options_copy[AZ_max_iter] = N_iterations;

   context = (struct aztec_context *) ML_allocate(sizeof(struct 
                                                         aztec_context));
   context->options        = options_copy;
   context->params         = params_copy;
   context->proc_config    = proc_config;
   context->comm           = ml_handle->comm;
   context->status         = status;
   context->prec_or_solver = N_iterations;
   matrix_name++;
   
   ML_Amat = &(ml_handle->Amat[level]);

   if (ML_Amat->matvec->ML_id == ML_EMPTY) {
       if (ml_handle->comm->ML_mypid == 0) {
          printf("ML_Gen_Smoother_Aztec: matvec not defined? \n");
       }
       exit(1);
   }

   if (global_mat_flag) {
      /* ----------------------------------------------------------------- */
      /* extract local matrix using getrow function and store it into a    */
      /* CSR data object                                                   */
      /* ----------------------------------------------------------------- */

      if ( level < 0 || level >= ml_handle->ML_num_levels ) {
         printf("ML_Gen_SmootherAztec error : invalid level number.\n");
         exit(-1);
      }
      op      = (ML_Operator *) &ml_handle->Amat[level];
      data    = op->data;
      osize   = op->outvec_leng;
      row_ptr = (int *) malloc(sizeof(int)*(osize+1));
      space   = osize * 5 + 30;
      getrow_flag = 0;
      if ( op->getrow->internal != NULL ) {
         getrow_flag = 1;
      } else if ( op->getrow->external != NULL ) {
         getrow_flag = 2;
      } else {
         printf("ML_Gen_CoarseSolverSuperLU error : no getrow function.\n");
         exit(-1);
      }
   
      flag    = 0;
   
      while (flag == 0) {
         cols    = (int    *) malloc(sizeof(int)*space);
         vals    = (double *) malloc(sizeof(double)*space);
   
         nz_ptr = 0;
         row_ptr[0] = nz_ptr;
         flag = 1;
         for (i = 0; i < osize; i++) {
            if ( getrow_flag == 1 ) {
               flag = op->getrow->internal((void*)op, 1, &i, space-nz_ptr,
                                 &(cols[nz_ptr]), &(vals[nz_ptr]), &length);
            } else {
               flag = op->getrow->external(data, 1, &i, space-nz_ptr,
                                  &(cols[nz_ptr]), &(vals[nz_ptr]), &length);
            }
   
            if (flag == 0) break;
            zero_flag = 1;
            for (j = 0; j < length; j++)
               if ( vals[nz_ptr+j] != 0.0 ) {zero_flag = 0; break;}
   
            if ( zero_flag == 1 )
            {
               cols[nz_ptr] = i;
               vals[nz_ptr] = 1.0;
               length = 1;
            }
            nz_ptr += length;
            row_ptr[i+1] = nz_ptr;
         }
         if (flag == 0) {
            dsize = (double) osize;
            di    = (double) (i+1);
            dsize = 1.2*dsize/di;
            space = (int) ( ((double) space)*dsize);
            space++;
            ML_free(vals);
            ML_free(cols);
         }
      }
      csr_mat = (ML_Matrix_DCSR *) malloc(sizeof(ML_Matrix_DCSR));
      csr_mat->mat_n  = osize;
      csr_mat->mat_ja = cols;
      csr_mat->mat_a  = vals;
      csr_mat->mat_ia = row_ptr;
      csr_mat->comminfo = op->getrow->pre_comm;
   
      /* ----------------------------------------------------------------- */
      /* form a global matrix                                              */
      /* ----------------------------------------------------------------- */
   
      csr2_mat = (ML_Matrix_DCSR *) malloc(sizeof(ML_Matrix_DCSR));
      ML_Gen_Amatrix_Global( csr_mat, csr2_mat, ml_handle->comm, &offset);
      free(row_ptr);
      free(cols);
      free(vals);
      free(csr_mat);
      nrows = csr2_mat->mat_n;
      AZ_Amat = AZ_matrix_create(nrows);
      sub_proc_config = (int *) malloc(AZ_PROC_SIZE*sizeof(int));
      sub_proc_config[AZ_node] = 0;
      sub_proc_config[AZ_N_procs] = 1;
   #ifdef ML_MPI
      tptr = AZ_get_comm(proc_config);
      AZ_set_comm(sub_proc_config, *tptr);
   #endif
      context->proc_config    = sub_proc_config;
      context->offset         = offset;
   
      AZ_set_MATFREE(AZ_Amat, csr2_mat, wrapper_DCSR_matvec);
      AZ_set_MATFREE_getrow(AZ_Amat, csr2_mat, wrapper_DCSR_getrow, NULL,
                     0, sub_proc_config);
      context->Amat  = AZ_Amat;
      AZ_Amat->must_free_data_org = 0;   /* data_org will be freed later by  */
                                         /* ML (via AZ_ML_SmootherClean) and */
                                         /* so Aztec does not need to free it*/
   }
   else {
   
      AZ_Amat = AZ_matrix_create(ML_Amat->invec_leng);
      context->Amat  = AZ_Amat;
      AZ_mlcomm2data_org(ML_Amat->getrow->pre_comm,&data_org);
      data_org[AZ_name] = matrix_name;
   
      if ((ML_Amat->matvec->ML_id == ML_EXTERNAL) &&
          (ML_Amat->matvec->external == az_matvec_wrapper)) {
   
         /* This matrix was originally generated by Aztec. A new   */ 
         /* data_org was made ... so that we could give the matrix */
         /* a new name to keep Aztec from getting confused.        */
   
         orig_context = (struct aztec_context *) ML_Amat->data;
         orig_Amat     = orig_context->Amat;
         orig_data_org = orig_Amat->data_org;
         data_org[AZ_matrix_type] = orig_data_org[AZ_matrix_type];
         data_org[AZ_N_internal]  = orig_data_org[AZ_N_internal];
         data_org[AZ_N_border]    = orig_data_org[AZ_N_border];
         data_org[AZ_N_int_blk]   = orig_data_org[AZ_N_int_blk];
         data_org[AZ_N_bord_blk]  = orig_data_org[AZ_N_bord_blk];
         data_org[AZ_N_ext_blk]   = orig_data_org[AZ_N_ext_blk];
   
         if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) 
            AZ_set_MSR(AZ_Amat, orig_Amat->bindx, orig_Amat->val, data_org,
                       0, NULL, AZ_LOCAL);
         else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) 
            AZ_set_VBR(AZ_Amat, orig_Amat->rpntr, orig_Amat->cpntr,
                       orig_Amat->bpntr, orig_Amat->indx,orig_Amat->bindx,
                       orig_Amat->val, data_org, 0, NULL, AZ_LOCAL);
         else {
            if (ml_handle->comm->ML_mypid == 0) {
               printf("AZ_set_ML_preconditioner: Can not use with");
               printf("Aztec matrix-free matrices\n");
            }
            exit(1);
         }
      }
      else if ((ML_Amat->matvec->ML_id == ML_INTERNAL) &&
               (ML_Amat->matvec->internal == MSR_matvec))  {
   
         /* This matrix was generated by ML  */ 
   
         data_org[AZ_matrix_type] = AZ_MSR_MATRIX;
         data_org[AZ_N_internal]  = 0;
         data_org[AZ_N_border]    = ML_Amat->invec_leng;
         data_org[AZ_N_int_blk]   = 0;
         data_org[AZ_N_bord_blk]  = ML_Amat->invec_leng;
         data_org[AZ_N_ext_blk]   = data_org[AZ_N_external];
   
         msr_mat        = (struct ML_CSR_MSRdata *) ML_Amat->data;
         AZ_set_MSR(AZ_Amat, msr_mat->columns, msr_mat->values, data_org,
                    0, NULL, AZ_LOCAL);
      }
      else {
         if (ml_handle->comm->ML_mypid == 0) {
             printf("Currently ML_Gen_Smoother_Aztec can only work with\n");
             printf("matrices generated by Aztec or ML's RAP feature\n");
         }
         exit(1);
      }
   }

   if (options_copy[AZ_precond] == AZ_user_precond) {
       if (prec_function == AZ_precondition) {
          if (ml_handle->comm->ML_mypid == 0) {
             printf("ML_Gen_SmootherAztec:");
             printf("\toptions[AZ_precond]=AZ_user_precond but \n");
             printf("\tprec_function is set to AZ_precondition (causing Aztec to\n");
             printf("\trecursively call AZ_precondition). Either set\n");
             printf("\toptions[AZ_precond] to another preconditioner\n");
             printf("\tor provide an alternative prec_function.\n");
          }
          exit(1);
       }
       if (prec_function == NULL) {
          if (ml_handle->comm->ML_mypid == 0) {
             printf("ML_Gen_SmootherAztec:");
             printf("\toptions[AZ_precond]=AZ_user_precond but \n");
             printf("\tprec_function is set to NULL(causing Aztec to\n");
             printf("\trecursively call AZ_precondition). Either set\n");
             printf("\toptions[AZ_precond] to another preconditioner\n");
             printf("\tor provide an alternative prec_function.\n");
          }
          exit(1);
       }
#ifdef AZ_ver2_1_0_3
       AZ_Prec = AZ_precond_create(context->Amat, prec_function, NULL);
#else
       AZ_Prec = AZ_precond_create(context->Amat, prec_function);
#endif
    }
#ifdef AZ_ver2_1_0_3
    else AZ_Prec = AZ_precond_create(context->Amat, AZ_precondition, NULL);
#else
    else AZ_Prec = AZ_precond_create(context->Amat, AZ_precondition);
#endif
    context->Prec = AZ_Prec;
#ifdef AZ_ver2_1_0_5
    context->scaling = AZ_scaling_create();
#endif
    ML_Set_Smoother(ml_handle,level,pre_or_post,(void *)context, 
                    az_wrap_solvers, "Aztec");

    /* hack in a function that will be invoked */
    /* by ML_Destroy() to clean up memory       */

    if (pre_or_post == ML_PRESMOOTH) 
       ml_handle->pre_smoother[level].data_destroy = AZ_ML_SmootherClean;
    else
       ml_handle->post_smoother[level].data_destroy= AZ_ML_SmootherClean;

   /* To use Aztec's preconditioners without using Aztec's solvers, */
   /* AZ_initialize must be called!!!                               */

   if (context->prec_or_solver == AZ_ONLY_PRECONDITIONER) {
      if (!AZ_initialize(NULL, NULL, context->options, context->params, 
		    context->status, context->proc_config, 
		    context->Amat, context->Prec, save_old_values
#ifdef AZ_ver2_1_0_5
                    ,context->scaling         
#endif
                    )) {
         exit(-1);
      }
#ifdef ML_EXPERIMENT
size = data_org[AZ_N_internal] + data_org[AZ_N_border];
xxx  = (double *) malloc( sizeof(double)*size );
rhss = (double *) malloc( sizeof(double)*size );
az_wrap_solvers(context, 1, xxx, 1, rhss);
free(rhss);
free(xxx);
#endif
   }

}

/*****************************************************************************/
/*****************************************************************************/

int az_wrap_solvers(void *data, int in, double x[], int out, 
                    double rhs[])
{
   struct aztec_context *context;
   int    *data_org, i, n, n2, one = 1;
   double *p2, alpha = 1.; 
   double temp, *global_rhs, *global_x, *orig_x;

   context = (struct aztec_context *) data;
   data_org = context->Amat->data_org;

   n        = data_org[AZ_N_internal] + data_org[AZ_N_border];
   n2       = n + data_org[AZ_N_external];
   p2       = (double *) AZ_allocate( (n2+1)*sizeof(double));

   if (p2 == NULL) {
      printf("az_wrap_solvers: Out of space\n"); exit(1);
   }

   /* we have replicated entire linear system on each processor */
   /* must gather rhs and initial guess                         */
   if ( n != in ) {
      ML_memory_alloc((void**) &global_rhs, n*sizeof(double),"LU1" );
      ML_memory_alloc((void**) &global_x,   n*sizeof(double),"LU2" );
      for ( i = 0; i <  n; i++ ) global_rhs[i] = 0.;
      for ( i = 0; i <  n; i++ ) global_x[i]   = 0.;
      for ( i = 0; i < in; i++ ) global_rhs[i] = rhs[i];
      for ( i = 0; i < in; i++ ) global_x[i]   = x[i];
      i = in; ML_Comm_GappendDouble(context->comm, global_rhs, &i, n);
      i = in; ML_Comm_GappendDouble(context->comm, global_x, &i, n);
      orig_x = x;
      x = global_x;
      rhs = global_rhs;
   }


   for (i = 0; i < n; i++) p2[i] = x[i];
   if (context->prec_or_solver == AZ_ONLY_PRECONDITIONER) {
      context->Amat->matvec(p2,x,context->Amat,context->proc_config);

      for (i = 0; i < n; i++) {
         temp  = p2[i];
         p2[i] = rhs[i] - x[i];
         x[i]  = temp;
      }
      context->Prec->prec_function(p2,context->options,
                                    context->proc_config,context->params,
                                    context->Amat, context->Prec);
      daxpy_(&n,&alpha, p2, &one, x, &one);
   }
   else {
      AZ_oldsolve(p2,rhs,context->options,context->params, 
                  context->status,context->proc_config,context->Amat,
                  context->Prec, context->scaling);
      for (i = 0; i < n; i++) x[i] = p2[i];
   }
   AZ_free(p2);

   /* we have replicated entire linear system on each processor */
   /* must gather solution                                      */
   if ( n != in ) {
      x = orig_x;
      for ( i = 0; i < in; i++ ) x[i] = global_x[i+context->offset];
      ML_memory_free((void**) &global_rhs);
      ML_memory_free((void**) &global_x );
   }
   return(1);
}

/*****************************************************************************/
/*****************************************************************************/

void AZ_ML_SmootherClean(void *data)
{
/*
 * Clean up space created by ML_Gen_SmootherAztec
 *
 */
   struct aztec_context *context;
   ML_Matrix_DCSR *csr2_mat;

   context = (struct aztec_context *) data;
   context->options[AZ_keep_info] = 0;

   /* we are using an Aztec smoother with an ML preconditoner */
   /* We must clean up the ml stuff.                          */
   if ( (context->options[AZ_precond] == AZ_user_precond) && 
        (context->Prec->prec_function == ML_precondition) &&
        (context->Prec->precond_data  != NULL) ) {
	ML_Solve_SmootherDestroy( context->Prec->precond_data );
   }




   AZ_iterate_finish(context->options, context->Amat, context->Prec);
   AZ_free(context->options); 
   AZ_free(context->params);
   AZ_free(context->Amat->data_org); 

   /* we have created a global matrix and must clean up storage */
   /* associated with this.                                     */
   /* we also have created a proc_config that must be cleaned.  */

   if (context->Amat->matvec == wrapper_DCSR_matvec) {
      csr2_mat = (ML_Matrix_DCSR *) context->Amat->matvec_data;
      ML_memory_free( (void **) &(csr2_mat->mat_ja));
      ML_memory_free( (void **) &(csr2_mat->mat_a) );
      ML_memory_free( (void **) &(csr2_mat->mat_ia) );
      free(csr2_mat);
      free(context->proc_config);
   }
   AZ_matrix_destroy(&(context->Amat) );
   AZ_precond_destroy(&(context->Prec));
#ifdef AZ_ver2_1_0_5
   AZ_scaling_destroy(&(context->scaling));
#endif
   AZ_free(context);
}
#endif

#ifdef AZTEC2_0

/*****************************************************************************/
/*****************************************************************************/

AZ_PRECOND *AZ_precond_create(AZ_MATRIX *Pmat, void (*prec_fun)(
        double *, int *, int *, double *, struct AZ_MATRIX_STRUCT  *,
               struct AZ_PREC_STRUCT *) )
{
   AZ_PRECOND *precond;

   precond = (AZ_PRECOND *) AZ_allocate(sizeof(AZ_PRECOND));
   if (precond == NULL) {
      printf("Error: Not enough space in AZ_precond_create().\n");
      exit(1);
   }
   precond->Pmat = Pmat;
   precond->prec_function = prec_fun;
   precond->options = NULL;
   precond->params  = NULL;
   return(precond);
}
AZ_MATRIX *AZ_matrix_create(int local)
{
/*
 * Create an Aztec AZ_MATRIX structure and fill in the noncommunication
 * related fields of data_org[].
 * Note: This matrix will not work properly with Aztec's AZ_exchange_bdry()
 *       subroutine. Instead, it is intended that this matrix be used for
 *       matrix-free users and matrices which do not require communication.
 *
 * Parameters
 * ========
 *   local              Number of matrix equations residing on this processor.
 *
 *   additional         local+additional is the required size of a vector, x,
 *                      that will be applied to a matrix when performing a
 *                      matrix-vector product. The first 'local' components of
 *                      'x' must contain the appropriate data. The remaining
 *                      'additional' components are used as workspace inside
 *                      the user's matrix vector product.
 *  matrix_type         Either AZ_MSR_MATRIX, AZ_VBR_MATRIX, or AZ_USER_MATRIX.
 *  local_blks          When matrix_type == AZ_VBR_MATRIX, 'local_blks'
 *                      indicates how many block equations reside on this node.
 */

   AZ_MATRIX  *Amat;

   Amat     = (AZ_MATRIX *) AZ_allocate(sizeof(AZ_MATRIX));
   if (Amat == NULL) {
      printf("Error: Not enough space in AZ_matrix_create().\n");
      exit(1);
   }
   Amat->matrix_type = AZ_none;
   Amat->rpntr       = NULL;
   Amat->cpntr       = NULL;
   Amat->bpntr       = NULL;
   Amat->bindx       = NULL;
   Amat->indx        = NULL;
   Amat->val         = NULL;
   Amat->data_org    = NULL;
   Amat->matvec      = NULL;
   Amat->matrix_norm = -1.0;
   Amat->aux_ival    = NULL;
   Amat->aux_dval    = NULL;
   Amat->aux_ptr     = NULL;
   Amat->aux_matrix  = NULL;

   return(Amat);
}

/*****************************************************************************/
/*****************************************************************************/

void AZ_set_MSR(AZ_MATRIX *Amat, int bindx[], double val[], int data_org[],
        int N_update, int update[], int option)
{
   Amat->bindx    = bindx;
   Amat->val      = val;
   Amat->data_org = data_org;
   Amat->matrix_type = AZ_MSR_MATRIX;
   Amat->matvec   = AZ_MSR_matvec_mult;
}

/*****************************************************************************/
/*****************************************************************************/

void AZ_set_VBR(AZ_MATRIX *Amat, int rpntr[], int cpntr[], int bpntr[],
        int indx[], int bindx[], double val[], int data_org[],
        int N_update, int update[], int option)
{
   Amat->rpntr = rpntr;
   Amat->cpntr = cpntr;
   Amat->bpntr = bpntr;
   Amat->indx  = indx;
   Amat->bindx = bindx;
   Amat->val   = val;
   Amat->data_org = data_org;
   Amat->matrix_type = AZ_VBR_MATRIX;
   Amat->matvec   = AZ_VBR_matvec_mult;
}

/*****************************************************************************/
/*****************************************************************************/

void AZ_matrix_destroy(AZ_MATRIX **Amat)
{
   AZ_free(*Amat);
   *Amat = NULL;
}

/*****************************************************************************/
/*****************************************************************************/

void AZ_precond_destroy(AZ_PRECOND **precond)
{
   AZ_free(*precond);
   *precond = NULL;
}



/*****************************************************************************/
/*****************************************************************************/

void AZ_set_proc_config(int proc_config[], MPI_AZComm comm)

{
  get_parallel_info(&(proc_config[AZ_node]), &(proc_config[AZ_N_procs]),
                    &(proc_config[AZ_dim]));
}


#endif

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

void AZ_mlcomm2data_org(ML_CommInfoOP *comm_info, int *data_org[])
{
/* NOTE: NOT EVERY comm_info op can be turned into a data_org */
   int i, j, count, *neighbors, *itemp, N_neighbors, total_send;
   int count2, *start_rcv = NULL, length, flag;

   N_neighbors = ML_CommInfoOP_Get_Nneighbors(comm_info);
   neighbors   = ML_CommInfoOP_Get_neighbors(comm_info);
   total_send = 0;
   if ( N_neighbors > AZ_MAX_NEIGHBORS) {
      printf("Need to increase AZ_MAX_NEIGHBORS in az_aztec_defs.h and \n");
      printf("recompile Aztec\n");
   }
   for (i = 0; i < N_neighbors; i++) {
      itemp  = ML_CommInfoOP_Get_rcvlist(comm_info, neighbors[i]);
      length = ML_CommInfoOP_Get_Nrcvlist(comm_info,neighbors[i]);
      if (itemp != NULL) {
          if (start_rcv == NULL) {
             start_rcv = ML_allocate((N_neighbors+1)*sizeof(int));
             if (start_rcv==NULL) pr_error("No space in AZ_mlcomm2data_org\n");
             for (j = 0; j < N_neighbors; j++) start_rcv[j] = -1;
          }
          /* check that receive list is contiguous (needed by Aztec) */
          flag = 0;
          for (j = 0; j < length-1; j++) 
             if ( itemp[j] != itemp[j+1]-1) flag = 1;
          if (flag == 1) {
             printf("AZ_mlcomm2data_org:I don't believe this comm object\n");
             printf("\t\twas created from RAP or Aztec\n");
             exit(1);
          }
          start_rcv[i] = itemp[0];
          free(itemp);
      }
      total_send += ML_CommInfoOP_Get_Nsendlist(comm_info,neighbors[i]);
   }
   if (start_rcv != NULL) {
      AZ_sort(start_rcv,N_neighbors, neighbors, NULL);
      ML_free(start_rcv);
   }

   *data_org = (int *) ML_allocate(((unsigned) total_send + AZ_send_list)
                                   *sizeof(int));
   if (*data_org == NULL) {
      (void) fprintf(stderr, "ERROR: Not enough dynamic space.\n");
      exit(-1);
   }
   (*data_org)[AZ_total_send] = total_send;

   count = AZ_send_list;
   count2 = 0;
   (*data_org)[AZ_N_neigh] = N_neighbors;
   for (i = 0; i < (*data_org)[AZ_N_neigh]; i++) {
        (*data_org)[AZ_neighbors+i] = neighbors[i];
        (*data_org)[AZ_send_length+i] = 
                    ML_CommInfoOP_Get_Nsendlist(comm_info,neighbors[i]);
        (*data_org)[AZ_rec_length+i] = 
                    ML_CommInfoOP_Get_Nrcvlist(comm_info,neighbors[i]);
        itemp = ML_CommInfoOP_Get_sendlist(comm_info, neighbors[i]);
        for (j = 0; j < (*data_org)[AZ_send_length+i]; j++)
           (*data_org)[count++] = itemp[j];
        free(itemp);
        count2 += (*data_org)[AZ_rec_length+i];
   }
   (*data_org)[AZ_N_external] = count2;
   free(neighbors);
}

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

void notusedAZML_convert_data_org(ML_Operator *matrix, int data_org[],
        int rcv_list[], int remap[], int leng, int add_or_not)
{
   int i, count, count2;



    ML_CommInfoOP_Set_neighbors( &(matrix->getrow->pre_comm),
                                 data_org[AZ_N_neigh],&(data_org[AZ_neighbors]),
                                 add_or_not, remap, leng);

    count = AZ_send_list;
    count2 = 0;

    if (rcv_list == NULL) {
       for (i = 0; i < data_org[AZ_N_neigh]; i++) {
          ML_CommInfoOP_Set_exch_info(matrix->getrow->pre_comm,
                    data_org[AZ_neighbors+i], data_org[AZ_rec_length+i], NULL,
                    data_org[AZ_send_length+i], &(data_org[count]));
          count += data_org[AZ_send_length+i];
       }
    }
    else {
       for (i = 0; i < data_org[AZ_N_neigh]; i++) {
          ML_CommInfoOP_Set_exch_info(matrix->getrow->pre_comm,
                    data_org[AZ_neighbors+i], data_org[AZ_rec_length+i],
                    &(rcv_list[count2]), data_org[AZ_send_length+i],
                    &(data_org[count]));
          count2 += data_org[AZ_rec_length+i];
          count += data_org[AZ_send_length+i];
       }
    }
}
int  wrapper_DCSR_getrow(int columns[], double values[], int row_lengths[],
        struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows,
        int requested_rows[], int allocated_space){

   int       status;
   ML_Matrix_DCSR *csr2_mat;
   struct ML_CSR_MSRdata *temp_ptr;

   csr2_mat = (ML_Matrix_DCSR *) Amat->getrow_data;

   temp_ptr = (struct ML_CSR_MSRdata *) malloc(sizeof(struct ML_CSR_MSRdata));
   temp_ptr->rowptr = csr2_mat->mat_ia;
   temp_ptr->columns= csr2_mat->mat_ja;
   temp_ptr->values = csr2_mat->mat_a;

   status = CSR_getrows(temp_ptr, N_requested_rows, requested_rows,
                        allocated_space, columns, values, row_lengths);

   free(temp_ptr);
   return(status);

}
void wrapper_DCSR_matvec(double *b, double *c,AZ_MATRIX *Amat,int proc_config[])
{
   ML_Matrix_DCSR *csr2_mat;
   struct ML_CSR_MSRdata *temp_ptr;

   if (proc_config[AZ_N_procs] > 1) {
      printf("wrapper_DCSR_getrow: only implemented in serial\n"); exit(1);
   }
   csr2_mat = (ML_Matrix_DCSR *) Amat->getrow_data;

   temp_ptr = (struct ML_CSR_MSRdata *) malloc(sizeof(struct ML_CSR_MSRdata));
   temp_ptr->rowptr = csr2_mat->mat_ia;
   temp_ptr->columns= csr2_mat->mat_ja;
   temp_ptr->values = csr2_mat->mat_a;

   localCSR_matvec(temp_ptr, csr2_mat->mat_n, b, csr2_mat->mat_n, c);
   free(temp_ptr);
}


#else

/* to satisfy the requirement of certain compilers */
int ML_empty;

#endif

