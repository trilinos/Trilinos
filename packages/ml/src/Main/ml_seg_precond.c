/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/************************************************************************/
/* Methods for segregated preconditioners                               */
/* Author: Dawn Chamberlain 12/3/99                                     */
/************************************************************************/

#ifdef ML_SEG
#include "ml_include.h"
#include "ml_aztec_utils.h"
#include "ml_seg_precond.h"

/************************************************************************/
/* Creates a SEG object that creates and stores submatrices of a        */
/* corresponding to an nblocks x nblocks block partition of A.  It is   */
/* intended to be used in a preconditioning scheme, so the SEG object   */
/* will only store the the diagonal, the upper triangle, or the lower   */
/* triangle (depending on the value of format)                          */
/*----------------------------------------------------------------------*/

int ML_SEG_Create(struct ML_SEG_Struct **seg, AZ_MATRIX *Amat, int nblocks, 
               int **rowlists, int *rowlengs, int format, int proc_config[])
{
   int i, j, noffdiags, nrows, ctr;

   (*seg) = (struct ML_SEG_Struct *) ML_allocate(sizeof(struct ML_SEG_Struct));
	
   (*seg)->SEG_nblocks = nblocks;
   (*seg)->SEG_format  = format;

   (*seg)->SEG_diag_list = 
       (ML_SEG_MATRIX **)ML_allocate(nblocks*sizeof(struct ML_SEG_MATRIX_Struct *));
   (*seg)->SEG_rowlist_lengs=(int *)ML_allocate(nblocks*sizeof(int));
   (*seg)->SEG_rowlists=(int **)ML_allocate(nblocks*sizeof(int *));

   if ((*seg)->SEG_rowlists == NULL) 
   {
      printf("ML_SEG_Create: memory allocation error\n");
      exit(-1);
   }
	
   /* allocate space for the rowlists and copy the values passed in */

   nrows=0;
   for (i = 0; i < nblocks; i++) 
   {
      (*seg)->SEG_rowlists[i]=(int *)ML_allocate(rowlengs[i]*sizeof(int));
      (*seg)->SEG_rowlist_lengs[i] = rowlengs[i];
      if ((*seg)->SEG_rowlists[i] == NULL) {
         printf("ML_SEG_Create: memory allocation error\n");
         exit(-1);
      }
		
      for (j = 0; j < rowlengs[i]; j++)
         (*seg)->SEG_rowlists[i][j]=rowlists[i][j];
	
      /* all formats require the diagonal, so create and store the */
      /* diagonal blocks: */

      (*seg)->SEG_diag_list[i] = ML_SEG_Matrix_Create(Amat, rowlengs[i], 
                                     rowlists[i], rowlengs[i], rowlists[i], 
                                     ML_SEG_DIAG_ELT, proc_config);
      nrows += rowlengs[i];
   }
		
   (*seg)->SEG_total_nrows = nrows;

   /* set the off diagonal list to null since we only care about the diagonal */

   if (format == ML_SEG_DIAGONAL) {
      (*seg)->SEG_offdiag_list = NULL;
      noffdiags = 0;
   }
   /* figure out how many blocks there are in the upper triangle and */
   /* allocate space */
   else if (format == ML_SEG_UPPER_TRIANGULAR) 
   {
      noffdiags = nblocks * (nblocks-1) / 2;
      (*seg)->SEG_offdiag_list = 
                  (ML_SEG_MATRIX **) ML_allocate(noffdiags*sizeof(ML_SEG_MATRIX *));
      ctr=0;
      /* for an upper triangular matrix, we will store the off diagonal */
      /* matrices starting with the bottom one (since when we apply the */
      /* preconditioner, we will start with the bottom block row and    */
      /* work our way up)  The ordering is as follows:
                   d 3 4 5
                   x d 1 2
                   x x d 0
		   x x x d
         where 'x's correspond to blocks we don't care about and 'd's 
         correspond to diagonal blocks (which are stored elsewhere) */
          
      for (i = nblocks-2; i >= 0; i--)
         for (j = i+1; j < nblocks; j++) 
         {
            (*seg)->SEG_offdiag_list[ctr] = 
                 ML_SEG_Matrix_Create(Amat, rowlengs[i], rowlists[i],
                   rowlengs[j], rowlists[j], ML_SEG_OFFDIAG_ELT, proc_config);
            ctr ++;
         }
   }
   else if (format == ML_SEG_LOWER_TRIANGULAR) 
   {
      noffdiags = nblocks * (nblocks-1) / 2;
      (*seg)->SEG_offdiag_list = 
            (ML_SEG_MATRIX **) ML_allocate(noffdiags*sizeof(ML_SEG_MATRIX *));
      ctr=0;
      /* for a lower triangular matrix, we will store the off diagonal */
      /* matrices starting with the top one (since when we apply the   */
      /* preconditioner, we will start with the top block row and work */
      /* our way down)  The ordering is as follows:
                   d x x x
                   0 d x x
                   1 2 d x
		   3 4 5 d
         where 'x's correspond to blocks we don't care about and 'd's 
         correspond to diagonal blocks (which are stored elsewhere) */

      for (i = 1; i < nblocks; i++)
         for (j = 0; j < i; j++) {
            (*seg)->SEG_offdiag_list[ctr] =   
               ML_SEG_Matrix_Create(Amat, rowlengs[i], rowlists[i],
                   rowlengs[j], rowlists[j], ML_SEG_OFFDIAG_ELT, proc_config);
            ctr++;
         }
   }
   else 
   {
      printf("ML_SEG_Create ERROR: Unknown format (%d)\n", format);
      exit(-1);
   }
   (*seg)->SEG_noffdiags = noffdiags;
   return(0);
} /* ML_SEG_Create */

/************************************************************************/
/* destructor                                                           */
/*----------------------------------------------------------------------*/

void ML_SEG_Destroy(struct ML_SEG_Struct **seg)
{
   int nblocks, noffdiags, i, **rowlists, *rowlistlengs;
   struct ML_SEG_MATRIX_Struct **diags, **offdiags;
	
   nblocks = (*seg)->SEG_nblocks;
   noffdiags = (*seg)->SEG_noffdiags;
   rowlists = (*seg)->SEG_rowlists;
   rowlistlengs = (*seg)->SEG_rowlist_lengs;
   diags = (*seg)->SEG_diag_list;
   offdiags = (*seg)->SEG_offdiag_list;

   for (i=0; i < nblocks; i++) ML_SEG_Matrix_Destroy(&(diags[i]));

   for (i=0; i < noffdiags; i++) ML_SEG_Matrix_Destroy(&(offdiags[i]));

   for (i=0; i < nblocks; i++) ML_free(rowlists[i]);

   ML_free(rowlistlengs);
   ML_free(rowlists);
   ML_free(diags);
   if (offdiags != NULL) ML_free(offdiags);
	
   ML_free(*seg);
} /* ML_SEG_Destroy */

/************************************************************************/
/* ML_SEG objects require additional information beyond a simple        */
/* AZ_MATRIX - a preconditioner can be associated with each SEG_MATRIX  */
/* (and hence each block of the original matrix can be preconditioned   */
/* separately).  When this function is called, submatrices are created, */
/* but if we later replace these submatrices, we will not copy those    */
/* matrices over (since that would be too much work).  Hence we need to */
/* keep track of whether or not each AZ matrix is the original one      */
/* created here - if we created it, we have to destroy it later, but if */
/* we didn't create it, we don't want to destroy it later               */
/*----------------------------------------------------------------------*/

struct ML_SEG_MATRIX_Struct *ML_SEG_Matrix_Create(AZ_MATRIX *Amat, int nrows, 
                                    int *rowlist, int ncols, int *collist, 
                                    int type, int proc_config[])
{
   struct ML_SEG_MATRIX_Struct *Smat;

   Smat = (struct ML_SEG_MATRIX_Struct *)
          ML_allocate(sizeof(struct ML_SEG_MATRIX_Struct));
   if (Smat == NULL) 
   {
      printf("ML_SEG_Matrix Create: allocation error\n");
      exit(-1);
   }
   Smat->AZ_options = NULL;
   Smat->AZ_params = NULL;
   Smat->AZ_scale = NULL;

   Smat->Amat = AZ_submatrix_create(Amat,nrows,rowlist,ncols,collist,proc_config);
   if (type == ML_SEG_DIAG_ELT) { /* only allow preconditioning on the diagonal */
      Smat->precond = AZ_precond_create(Smat->Amat, AZ_precondition, NULL);
   }
   else Smat->precond = NULL;

   Smat->Amat_changed = AZ_FALSE;
   Smat->precond_changed = AZ_FALSE;

   return (Smat);
} /*ML_SEG_Matrix_Create */

/************************************************************************/
/* matrix destructor                                                    */
/*----------------------------------------------------------------------*/

void ML_SEG_Matrix_Destroy(struct ML_SEG_MATRIX_Struct **Smat)
{
   AZ_PRECOND *precond = (*Smat)->precond;

   /* if this matrix hasn't been changed since we created it,
		 we have to destroy it */
   if ((*Smat)->Amat_changed == AZ_FALSE)
      AZ_submatrix_destroy(&((*Smat)->Amat));

   if ((precond != NULL) && ((*Smat)->precond_changed == AZ_FALSE)) 
   {
      AZ_free((void *) precond->context->tag);
      AZ_free((void *) precond->context);
      AZ_precond_destroy(&precond);
   }
	
   if ((*Smat)->AZ_options != NULL) ML_free((*Smat)->AZ_options);

   if ((*Smat)->AZ_params != NULL) ML_free((*Smat)->AZ_params);
	
   if((*Smat)->AZ_scale != NULL) AZ_scaling_destroy(&((*Smat)->AZ_scale));
	
   ML_free(*Smat);
} /* ML_SEG_Datrix_Destroy */


/************************************************************************/
/* set SEG preconditioner in Aztec                                      */
/*----------------------------------------------------------------------*/

void AZ_Set_SEG_Preconditioner(AZ_PRECOND *Precond,struct AZ_SEG_Struct *seg,
			      int options[])
{
   options[AZ_precond]    = AZ_user_precond;
   /* use the ml_ptr to store seg because it is convenient */
   Precond->ml_ptr = (void *) seg;
   Precond->prec_function = ML_SEG_Precondition;
} /* AZ_Set_SEG_Preconditioner */
 
/************************************************************************/
/* This is the analog of AZ_precondition - it calls the preconditioners */
/* for each diagonal block individually and multiplies in the           */
/* offdiagonal blocks (or the alternative matrices provided for the     */
/* offdiagonal blocks) as requested                                     */
/* Preconditioning wrapper function to be called by Aztec when using    */
/* a segregated solver as a preconditioner.                             */
/*----------------------------------------------------------------------*/

void ML_SEG_Precondition(double ff[], int options[], int proc_config[], 
                      double params[], AZ_MATRIX *mat, AZ_PRECOND *prec)
{
   int         i, j, k, lenf, nblocks, ctr, *rowlist, nrows, *tmp_opt;
   double      *fftmp, *vect, *tmp_par;
   struct ML_SEG_Struct *seg;
   struct ML_SEG_MATRIX_Struct *curMat, *offdiagMat;
   AZ_MATRIX  *Amat;
   AZ_PRECOND *precond;

   seg     = (struct ML_SEG_Struct *) prec->ml_ptr;
   nblocks = seg->SEG_nblocks;
   lenf    = seg->SEG_total_nrows;
	
   /* need temporary vectors for storing the vectors that go into and */
   /* come out of the calls to each block.  We know no single block   */
   /* will be larger than the whole matrix, so this is a convenient   */
   /* upper bound */

   vect = (double *)ML_allocate(lenf*sizeof(double));
   fftmp = (double *)ML_allocate(lenf*sizeof(double));
   if (fftmp == NULL) 
   {
      printf("ML_SEG_precondition: memory allocation error\n");
      exit(-1);
   }
   if (seg->SEG_format == ML_SEG_DIAGONAL) 
   {
      for (i=0; i < nblocks; i++) 
      {
         curMat = seg->SEG_diag_list[i];
         precond = curMat->precond;
         Amat = curMat->Amat;
         /* AZ_precond hasn't been set to AZ_none (0), precondition each */
         /* diagonal block.  If AZ_precond = AZ_none, assume the user    */
         /* doesn't want any preconditioning in this solve */
         if (options[AZ_precond]) 
         {    
            /* copy the relavent rows of the input vector into a vector we can
               send to the preconditioner */
            for (j=0; j < seg->SEG_rowlist_lengs[i]; j++) 
               fftmp[j] = ff[seg->SEG_rowlists[i][j]];
            /* if different options  have been set for this, use them.  Otherwise
               use the options sent into the overall preconditioner */
            if (curMat->AZ_options == NULL) tmp_opt = options;
            else tmp_opt = curMat->AZ_options;
            if (curMat->AZ_params == NULL) tmp_par = params;
            else tmp_par = curMat->AZ_params;
				
            precond->prec_function(fftmp,tmp_opt,proc_config,tmp_par,Amat,precond);

            /* now copy the preconditioned values into the vector that will */
            /* be returned */
            for (j=0; j < seg->SEG_rowlist_lengs[i]; j++) 
               ff[seg->SEG_rowlists[i][j]] = fftmp[j];
         }
      }
   }
   else if (seg->SEG_format == ML_SEG_UPPER_TRIANGULAR)
   {
      ctr = 0;
      for (i = nblocks-1; i >= 0; i--) {
         curMat = seg->SEG_diag_list[i];
         precond = curMat->precond;
         Amat = curMat->Amat;
         nrows = seg->SEG_rowlist_lengs[i];
         rowlist = seg->SEG_rowlists[i];
         /* deal with the diagonal as above */
         if (options[AZ_precond]) 
         {    
            for (j=0; j < seg->SEG_rowlist_lengs[i]; j++) 
            fftmp[j] = ff[seg->SEG_rowlists[i][j]];

            if (curMat->AZ_options == NULL) tmp_opt = options;
            else tmp_opt = curMat->AZ_options;

            if (curMat->AZ_params == NULL) tmp_par = params;
            else tmp_par = curMat->AZ_params;

            precond->prec_function(fftmp, tmp_opt, proc_config, tmp_par, 
                                   Amat, precond); 
            for (j=0; j < seg->SEG_rowlist_lengs[i]; j++) 
               ff[seg->SEG_rowlists[i][j]] = fftmp[j];
         }
         /* Now multiply in the offdiagonal blocks corresponding to this */
         /* block row.  We will go through them in order (since they     */
         /* were put there in the proper order).                         */ 

         for (j = i+1; j < nblocks; j++) 
         {
            offdiagMat = seg->SEG_offdiag_list[ctr++];
            /* this time we need to use 2 separate vectors for the data */
            /* going into and coming out of matvec since matvec does    */
            /* not overwrite its input the way the preconditioners do.  */
            for (k = 0; k < seg->SEG_rowlist_lengs[j]; k++)
               vect[k]=ff[seg->SEG_rowlists[j][k]];
            offdiagMat->Amat->matvec(vect, fftmp, offdiagMat->Amat, proc_config);
            for (k = 0; k < nrows; k++) ff[rowlist[k]] += fftmp[k];
         }
      }
   }
   else if (seg->SEG_format == ML_SEG_LOWER_TRIANGULAR) 
   {
      /* this section is completely analagous to the ML_SEG_UPPER_TRIANGULAR 
          case, except we go through the block matrix from top to 
          bottom instead of from bottom to top */
      ctr = 0;
      for (i = 0; i < nblocks ; i++) 
      {
         curMat = seg->SEG_diag_list[i];
         precond = curMat->precond;
         Amat = curMat->Amat;
         nrows = seg->SEG_rowlist_lengs[i];
         rowlist = seg->SEG_rowlists[i];
         if (options[AZ_precond]) 
         {    
            for (j=0; j < seg->SEG_rowlist_lengs[i]; j++) 
               fftmp[j] = ff[seg->SEG_rowlists[i][j]];
            if (curMat->AZ_options == NULL) tmp_opt = options;
            else tmp_opt = curMat->AZ_options;

            if (curMat->AZ_params == NULL) tmp_par = params;
            else tmp_par = curMat->AZ_params;

            precond->prec_function(fftmp, tmp_opt, proc_config, tmp_par, 
                                   Amat, precond);

            for (j=0; j < seg->SEG_rowlist_lengs[i]; j++) 
               ff[seg->SEG_rowlists[i][j]] = fftmp[j];
         }
         for (j = 0; j < i; j++) 
         {
            offdiagMat = seg->SEG_offdiag_list[ctr++];
            for (k = 0; k < seg->SEG_rowlist_lengs[j]; k++)
               vect[k]=ff[seg->SEG_rowlists[j][k]];
			
            offdiagMat->Amat->matvec(vect,fftmp,offdiagMat->Amat,proc_config);

            for (k = 0; k < nrows; k++) ff[rowlist[k]] += fftmp[k];
         }
      }
   }
   else 
   {
      printf("ML_SEG_precondition error:  Unknown matrix format (%d)\n", 
             seg->SEG_format);
      exit(-1);
   }
	
   ML_free(fftmp);
   ML_free(vect);
} /* ML_SEG_Precondition */


/************************************************************************/
/* This function sets an ml preconditioner for a block of a segregated  */
/* matrix.  The ml object must be defined in the calling function       */
/* (which must be generated using the appropriate submatrix block of A) */
/*----------------------------------------------------------------------*/

void ML_SEG_Set_ML_Precond(struct ML_SEG_Struct *seg, ML *ml, int block_row, 
	                   double params[], int options[], int proc_config[])
{
   /* note: ML preconditioners can only be specified for diagonal blocks */
   int *prec_options, i, save_old_values[6];
   double status[AZ_STATUS_SIZE];
   AZ_PRECOND *precond;
   struct ML_SEG_MATRIX_Struct *Smat;

   if (block_row > seg->SEG_nblocks) 
   {
      printf("ML_SEG_Set_ML_Precond error: requested block %d in a matrix \
             with %d blocks\n", block_row, seg->SEG_nblocks);
      exit(-1);
   }

   /* copy the options sent in and store them in seg->AZ_opt since aztec */
   /* expects certain options from an ml preconditioner */

   prec_options = (int *)ML_allocate(AZ_OPTIONS_SIZE * sizeof(int));
   if (prec_options == NULL) 
   {
      printf("memory allocation error in ML_SEG_Set_ML_Precond\n");
      exit(-1);
   }
	
   for (i = 0; i < AZ_OPTIONS_SIZE; i++) prec_options[i] = options[i];
	
   Smat = seg->SEG_diag_list[block_row];
   precond = Smat->precond; 
   AZ_precond_destroy( &(Smat->precond) );
   AZ_set_ML_preconditioner(&(Smat->precond), Smat->Amat, ml, prec_options);

   /* since we don't call this preconditioner directly from AZ_solve, */
   /* we have to initialize it ourselves */

   Smat->AZ_scale = AZ_scaling_create();

   if (!AZ_initialize(NULL, NULL, prec_options, params, status, proc_config, 
              Smat->Amat, Smat->precond, save_old_values, Smat->AZ_scale)) 
   {
      exit(-1);
   }

   Smat->AZ_options = prec_options;
} /* ML_SEG_Seg_ML_Precond */

/************************************************************************/
/* Converts a SEG matrix into an ML matrix (needed because the blocks   */
/* of SEG are based on AZ matrices)                                     */
/*----------------------------------------------------------------------*/
			
void ML_SEG_ML_Set_Amat(struct ML_SEG_Struct *seg, ML *ml, int level, 
                     int block_row, int proc_config[])
{
   int leng;
   struct ML_SEG_MATRIX_Struct *Smat;
   AZ_MATRIX *Amat;
	
   if (block_row > seg->SEG_nblocks) 
   {
      printf("ML_SEG_ML_SetAmat ERROR: requested block %d in a matrix \
             with %d blocks\n", block_row, seg->SEG_nblocks);
      exit(-1);
   }
   leng = seg->SEG_rowlist_lengs[block_row];
	
   Smat = seg->SEG_diag_list[block_row];
   Amat = Smat->Amat;

   AZ_ML_Set_Amat(ml, level, leng, leng, Amat, proc_config);
} /* ML_SEG_ML_Set_Amat */

/************************************************************************/
/* Function for replacing a submatrix of seg with another AZ_MATRIX.    */
/* Useful for instance in Wathen's scheme, where Mp is used instead of  */
/* the (2,2) block of A to precondition those variables                 */
/*----------------------------------------------------------------------*/

void ML_SEG_Replace_Submat(struct ML_SEG_Struct *seg, int block_row, 
                    int block_col, AZ_MATRIX *newMat, int proc_config[])
{
   int nblocks, block_loc;
   AZ_MATRIX *oldMat;
   AZ_PRECOND *oldPrec, *newPrec;
   struct ML_SEG_MATRIX_Struct *Smat;

   /* check if the matrix to be replaced is on the diagonal or not (and */
   /* if it is a valid block) and locate the block in the appropriate   */
   /* list in seg  */

   nblocks=seg->SEG_nblocks;
	
   if (block_row == block_col) 
      if (block_row < nblocks) Smat = seg->SEG_diag_list[block_row];
      else {
        printf("Can't replace diagonal block %d in a %d by %d block matrix\n", 
                block_row, nblocks, nblocks);
        exit(-1);
      }
   else if (seg->SEG_format == ML_SEG_UPPER_TRIANGULAR) 
   {
      if ((block_row >= nblocks-1) || (block_col <= block_row)) {
         printf("This is an upper triangular matrix - can't replace a block \
                in the lower triangle\n");
         printf("(block_row=%d and block_col=%d, nblocks=%d\n", block_row, 
                block_col, nblocks);
         exit(-1);
      }
      block_loc = (nblocks-block_row-1)*(nblocks-block_row-2)/2 + 
                  block_col-block_row-1;                                
      /* #blocks below this row and #blocks to left of this col */
      Smat = seg->SEG_offdiag_list[block_loc];
   }
   else if (seg->SEG_format == ML_SEG_LOWER_TRIANGULAR) 
   {
      nblocks=seg->SEG_nblocks;
      if ((block_col >= nblocks-1) || (block_col >= block_row)) 
      {
         printf("This is a lower triangular matrix - can't replace a block \
                in the upper triangle\n");
         printf("(block_row=%d and block_col=%d, nblocks=%d\n", block_row, 
                block_col, nblocks);
         exit(-1);
      }
      block_loc = block_row * (block_row-1)/2 /* #blocks above this row */
		  + block_col;   /* #blocks to left of this col */
      Smat = seg->SEG_offdiag_list[block_loc];
   }

   /* check if this block has already been changed, and free the old one if */
   /* it hasn't */
   if (Smat->Amat_changed == AZ_FALSE) 
   {
      oldMat = Smat->Amat;
      AZ_submatrix_destroy(&oldMat);
   }

   /* put this matrix in the appropriate place and mark that we have changed it */
   Smat->Amat = newMat;
   Smat->Amat_changed = AZ_TRUE;

   /* if this is a diagonal block, we also have to replace the preconditioner */
   /* associated with this matrix and free the old one */

   if (block_row == block_col) 
   {
      oldPrec = Smat->precond;
      newPrec = AZ_precond_create(newMat, AZ_precondition, NULL);
      AZ_precond_destroy(&oldPrec);
      Smat->precond = newPrec;
   }
} /* ML_SEG_Replace_Submat */

/************************************************************************/
/* An exact analog to AZ_Set_ML_Precond                                 */
/*----------------------------------------------------------------------*/

void ML_SEG_Set_AZ_Precond(struct ML_SEG_Struct *seg,int block_row,
                           double params[], int options[], int proc_config[])
{
   int *prec_opts, i, save_old_values[6];
   double status[AZ_STATUS_SIZE], *prec_pars;
   AZ_MATRIX *Amat;
   AZ_PRECOND *precond;
   struct ML_SEG_MATRIX_Struct *Smat;

   if (block_row > seg->SEG_nblocks) 
   {
      printf("ML_SEG_Set_AZ_Precond error: requested block %d in a matrix \
             with %d blocks\n", block_row, seg->SEG_nblocks);
      exit(-1);
   }
	
   prec_opts = (int *)ML_allocate(AZ_OPTIONS_SIZE * sizeof(int));
   prec_pars = (double *)ML_allocate(AZ_PARAMS_SIZE * sizeof(double));
   if (prec_pars == NULL) 
   {
      printf("memory alloction error in ML_SEG_Set_AZ_Precond\n");
      exit(-1);
   }

   for (i=0; i < AZ_OPTIONS_SIZE; i++) prec_opts[i] = options[i];
   for (i=0; i < AZ_PARAMS_SIZE; i++) prec_pars[i] = params[i];

   Smat = seg->SEG_diag_list[block_row];
   precond = Smat->precond;
   Amat = Smat->Amat;

   /*tmp = AZ_set_solver_parameters(params, prec_opts, Amat, precond); */
   Smat->AZ_scale = AZ_scaling_create();

   if (!AZ_initialize(NULL, NULL, prec_opts, params, status, proc_config,
                      Smat->Amat, precond, save_old_values, Smat->AZ_scale)) 
   {
      exit(-1);
   }
   Smat->AZ_options = prec_opts;
   Smat->AZ_params = prec_pars;
} /* ML_SEG_Set_AZ_Precond */

#else
/* to satisfy the requirement of certain compilers */
int ML_empty2;

#endif








