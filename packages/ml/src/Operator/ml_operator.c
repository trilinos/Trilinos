/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for the ML_Operator structure                              */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : March, 1999                                          */
/* ******************************************************************** */

#include "ml_operator.h"
#include <string.h>

/************************************************************************/
/* Create an ML_matrix and initialize relevant fields.                  */
/************************************************************************/

ML_Operator *ML_Operator_Create(ML_Comm *comm)
{
   ML_Operator *temp;

   temp = (ML_Operator *) ML_allocate(sizeof(ML_Operator));
   ML_Operator_Init(temp,comm);

   return(temp);
}

/************************************************************************/
/* destructor                                                           */
/************************************************************************/

int ML_Operator_Destroy( ML_Operator *mat)
{
   ML_Operator_Clean(mat);
   ML_free(mat);

   return 0;
}

/* ******************************************************************** */
/* Initialize                                                           */
/* ******************************************************************** */

int ML_Operator_Init( ML_Operator *mat, ML_Comm *comm)
{
   mat->ML_id = ML_ID_OP;
   ML_memory_alloc((void**)&(mat->matvec),sizeof(ML_Function),"OF1");
   mat->matvec->ML_id    = ML_EMPTY;
   mat->matvec->Nrows    = 0;
   mat->matvec->internal = NULL;
   mat->matvec->external = NULL;
   ML_memory_alloc((void**)&(mat->getrow),sizeof(ML_GetrowFunc),"OF2");
   mat->getrow->ML_id            = ML_EMPTY;
   mat->getrow->Nrows            = 0;
   mat->getrow->pre_comm         = NULL;
   mat->getrow->post_comm        = NULL;
   mat->getrow->internal         = NULL;
   mat->getrow->external         = NULL;
   mat->getrow->data             = NULL;
   mat->getrow->use_loc_glob_map = ML_NO;
   mat->getrow->loc_glob_map     = NULL;
   mat->getrow->row_map          = NULL;

   mat->to                  = NULL;
   mat->from                = NULL;
   mat->invec_leng          = 0;
   mat->outvec_leng         = 0;
   mat->data                = NULL;
   mat->diagonal            = NULL;      
   mat->N_nonzeros          = 0;
   mat->max_nz_per_row      = 0;
   mat->sub_matrix          = NULL;
   mat->from_an_ml_operator = 0;
   mat->data_destroy        = NULL;
   mat->build_time          = 0.0;
   mat->apply_time          = 0.0;
   mat->label               = NULL;
   mat->comm                = comm;
   mat->num_PDEs            = 1;
   return 0;
}

/* ******************************************************************** */
/* Clean (corresponding to Init                                         */
/* ******************************************************************** */

int ML_Operator_Clean( ML_Operator *mat)
{
#ifdef ML_DETAILED_TIMING
   double t1;

   if ( (mat->label != NULL) && ( mat->build_time != 0.0)) 
   {
      t1 = ML_gsum_double(mat->build_time, mat->comm);
      t1 = t1/((double) mat->comm->ML_nprocs);
      if ( (mat->comm->ML_mypid == 0) && (t1 != 0.0))
         printf(" Build time for %s (average) \t= %e\n",mat->label,t1);
      t1 = ML_gmax_double(mat->build_time, mat->comm);
      if ( (mat->comm->ML_mypid == 0) && (t1 != 0.0))
         printf(" Build time for %s (maximum) \t= %e\n",mat->label,t1);
      t1 = - mat->build_time;
      t1 = ML_gmax_double(t1, mat->comm);
      t1 = - t1;
      if ( (mat->comm->ML_mypid == 0) && (t1 != 0.0))
         printf(" Build time for %s (minimum) \t= %e\n",mat->label,t1);
   }
   if  (mat->label != NULL) {
      t1 = ML_gsum_double(mat->apply_time, mat->comm);
      t1 = t1/((double) mat->comm->ML_nprocs);
      if ( (mat->comm->ML_mypid == 0) && (t1 != 0.0))
         printf(" Apply time for %s (average) \t= %e\n",mat->label,t1);
      t1 = ML_gmax_double(mat->apply_time, mat->comm);
      if ( (mat->comm->ML_mypid == 0) && (t1 != 0.0))
         printf(" Apply time for %s (maximum) \t= %e\n",mat->label,t1);
      t1 = - mat->apply_time;
      t1 = ML_gmax_double(t1, mat->comm);
      t1 = - t1;
      if ( (mat->comm->ML_mypid == 0) && (t1 != 0.0))
         printf(" Apply time for %s (minimum) \t= %e\n",mat->label,t1);
   }
#endif

   if (mat->sub_matrix != NULL) ML_Operator_Destroy(mat->sub_matrix);
   if ((mat->data_destroy != NULL) && (mat->data != NULL)) {
		 /*printf("ready to call destroy %u\n",mat->data); */
      mat->data_destroy(mat->data);
      mat->data = NULL;
   }
   if (mat->diagonal != NULL) {
      ML_DVector_Destroy( &(mat->diagonal) );
   }

   mat->matvec->ML_id  = ML_ID_DESTROYED;
   mat->getrow->ML_id  = ML_ID_DESTROYED;
   if (mat->getrow->row_map != NULL) ML_free(mat->getrow->row_map);
   ML_CommInfoOP_Destroy(mat->getrow->pre_comm);
   ML_CommInfoOP_Destroy(mat->getrow->post_comm);

   if (mat->getrow->loc_glob_map != NULL) 
      ML_free(mat->getrow->loc_glob_map);
   ML_memory_free((void**)&(mat->matvec));
   ML_memory_free((void**)&(mat->getrow));
   if (mat->label != NULL) { free(mat->label); mat->label = NULL; }
   mat->num_PDEs            = 1;

   return 0;
}

/* ******************************************************************** */
/* The `half' functions are used to make a copy of a matrix without     */
/* copying the real data. This is used primarily when we want to square */
/* a matrix. 2matmult() makes changes to the input arguments            */
/* (temporarily). Thus, if we repeat the same argument twice, things    */
/* don't work.                                                          */
/* ******************************************************************** */

ML_Operator *ML_Operator_halfClone( ML_Operator *original)
{
   ML_Operator *mat;

   mat = ML_Operator_Create(original->comm);

   mat->ML_id = ML_ID_OP;
   mat->matvec->ML_id    = original->matvec->ML_id;
   mat->matvec->Nrows    = original->matvec->Nrows;
   mat->matvec->internal = original->matvec->internal;
   mat->matvec->external = original->matvec->external;
   mat->getrow->ML_id            = original->getrow->ML_id;
   mat->getrow->Nrows            = original->getrow->Nrows;
   mat->getrow->pre_comm         = original->getrow->pre_comm;
   mat->getrow->post_comm        = original->getrow->post_comm;
   mat->getrow->internal         = original->getrow->internal;
   mat->getrow->external         = original->getrow->external;
   mat->getrow->data             = original->getrow->data;
   mat->getrow->use_loc_glob_map = original->getrow->use_loc_glob_map;
   mat->getrow->loc_glob_map     = original->getrow->loc_glob_map;
   mat->getrow->row_map          = original->getrow->row_map;

   mat->to                  = original->to;
   mat->from                = original->from;
   mat->invec_leng          = original->invec_leng;
   mat->outvec_leng         = original->outvec_leng;
   mat->data                = original->data;
   mat->diagonal            = original->diagonal;
   mat->N_nonzeros          = original->N_nonzeros;
   mat->max_nz_per_row      = original->max_nz_per_row;
   mat->sub_matrix          = original->sub_matrix;
   mat->from_an_ml_operator = original->from_an_ml_operator;
   mat->data_destroy        = NULL;
   mat->build_time          = original->build_time;
   mat->apply_time          = original->apply_time;
   mat->label               = original->label;
   mat->comm                = original->comm;
   mat->num_PDEs            = original->num_PDEs;
   return mat;
}

/* ******************************************************************** */
/* destructor corresponding to halfClone                                */
/* ******************************************************************** */

int ML_Operator_halfDestroy( ML_Operator *mat)
{
   mat->sub_matrix = NULL;
   mat->diagonal   = NULL;
   mat->getrow->row_map = NULL;
   mat->getrow->loc_glob_map = NULL;
   mat->getrow->pre_comm = NULL;
   mat->getrow->post_comm = NULL;
   mat->label = NULL;
   return(ML_Operator_Destroy(mat));
}

/* ******************************************************************** */
/* Set the to and from field of the ML_Operator data structure          */
/* ******************************************************************** */

int ML_Operator_Set_1Levels(ML_Operator *mat,ML_1Level *from,ML_1Level *to)
{
   if ( mat->ML_id != ML_ID_OP ) {
      printf("ML_Operator_Set_1Levels error : wrong object.\n");
      exit(-1);
   }
   mat->to   = to;
   mat->from = from;
   return 0;
}

/* ******************************************************************** */
/* Set the BCs field of the ML_Operator data structure                  */
/* ******************************************************************** */

int ML_Operator_Set_BdryPts(ML_Operator *mat, ML_BdryPts *bc)
{
   if ( mat->ML_id != ML_ID_OP ) {
      printf("ML_Operator_Set_BdryPts error : wrong object.\n");
      exit(-1);
   }
   mat->bc = bc;
   return 0;
} 

/* ******************************************************************** */
/* Set the matvec information                                           */
/* ******************************************************************** */

int ML_Operator_Set_ApplyFuncData(ML_Operator *mat, int inlen, int outlen,
            int type, void *data, int nrows, 
            int (*func)(void*,int,double*,int,double*), int flag)
{
   if ( mat->ML_id != ML_ID_OP ) {
      printf("ML_Operator_Set_ApplyFunc error : wrong object.\n");
      exit(-1);
   }
/* newly added : 8/17/00 */
   if ( mat->data != NULL && mat->data_destroy != NULL ) 
   {
      mat->data_destroy(mat->data);
      mat->data = NULL;
   }
   mat->invec_leng = inlen;
   mat->outvec_leng = outlen;
   mat->data = data;
   if ( type == ML_INTERNAL ) mat->matvec->internal = func;
   else                       mat->matvec->external = func;
   mat->matvec->ML_id = type;
   mat->matvec->Nrows = nrows;
   if ( flag != 0 ) mat->from_an_ml_operator = flag;
   return 0;
}

/* ******************************************************************** */
/* Set the matvec information                                           */
/************************************************************************/

int ML_Operator_Set_ApplyFunc(ML_Operator *Op, int internal_or_external,
                       int (*func)(void *, int, double *, int, double *))
{
   if (internal_or_external == ML_EXTERNAL)
        Op->matvec->external = func;
   else Op->matvec->internal = func;

   Op->matvec->ML_id = internal_or_external;
   return 0;
}

/* ******************************************************************** */
/* Set matrix diagonal                                                  */
/************************************************************************/

int ML_Operator_Set_Diag(ML_Operator *Op, int size, double diagonal[])
{
   if (Op->diagonal != NULL) {
     pr_error("ML_Operator_Set_Diagonal: Diagonal is already nonNull. \
               It appears that diagonal already exists!\n");
   }
   ML_DVector_Create( &(Op->diagonal), NULL );
   ML_DVector_LoadData( Op->diagonal, size, diagonal );
   if (Op->outvec_leng != size)
     pr_error("ML_Operator_Set_Diagonal: Size (%d) does not match matrix \
        outvec length (%d)\n",size,Op->outvec_leng);
   return 0;
}

/* ******************************************************************** */
/* set getrow function                                                  */
/* ******************************************************************** */

int ML_Operator_Set_Getrow(ML_Operator *Op, int internal_or_external,
        int size, int (*func)(void *,int,int*,int,int*,double*,int*))
{
   if (internal_or_external == ML_EXTERNAL)
        Op->getrow->external = func;
   else Op->getrow->internal = func;

   Op->getrow->ML_id = internal_or_external;
   Op->getrow->Nrows = size;

   return 0;
}

/* ******************************************************************** */
/* get a requested row from the operator                                */
/* ******************************************************************** */

int ML_Operator_Getrow(ML_Operator *Amat, int N_requested_rows, 
                int requested_rows[], int allocated_space, int columns[], 
                double values[], int row_lengths[])
{
   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("ML_Operator_Getrow : Amat getrow not defined\n");

   if (Amat->getrow->ML_id == ML_EXTERNAL)
      return(Amat->getrow->external(Amat->data,N_requested_rows,
	requested_rows, allocated_space, columns, values, row_lengths));
   else 
      return(Amat->getrow->internal(Amat,N_requested_rows,requested_rows, 
                          allocated_space, columns, values, row_lengths));
}

/* ******************************************************************** */
/* get matrix diagonal                                                  */
/* ******************************************************************** */

int ML_Operator_Get_Diag(ML_Operator *Amat, int length, double **diag)
{
   int allocated_space, *cols, i, j, n;
   double *vals, *tdiag;

   if (Amat->diagonal == NULL)
   {
      if (Amat->getrow->ML_id == ML_EMPTY)
         pr_error("Error(ML_Operator_Get_Diag): diagonal not available\n");
      else
      {
         allocated_space = 30;
         cols = (int    *) malloc(allocated_space*sizeof(int   ));
         vals = (double *) malloc(allocated_space*sizeof(double));
         tdiag = (double *) malloc(length*sizeof(double));
         if (tdiag == NULL) 
            pr_error("Error(ML_Operator_Get_Diag): not enough space\n");
         for (i = 0; i < length; i++) tdiag[i] = 0.;
         for (i = 0; i < length; i++)
         {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,
                                     cols,vals,&n) == 0)
            {
               allocated_space = 2*allocated_space + 1;
               free(vals); free(cols);
               cols = (int    *) malloc(allocated_space*sizeof(int   ));
               vals = (double *) malloc(allocated_space*sizeof(double));
               if (vals == NULL)
               {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < n; j++)
               if (cols[j] == i) tdiag[i] = vals[j];
         }
         free(cols); free(vals);
         ML_Operator_Set_Diag(Amat, length, tdiag);
         free(tdiag);
      }
   }
   ML_DVector_GetDataPtr( Amat->diagonal, diag);
   return 0;
}


/* ******************************************************************** */
/* apply the operator to a vector                                       */
/************************************************************************/

int ML_Operator_Apply(ML_Operator *Op, int inlen, double din[], int olen,
                      double dout[])
{
#ifdef ML_TIMING
   double t0;

   t0 = GetClock();
#endif
   if (Op->matvec->ML_id == ML_EMPTY)
      pr_error("ML_Operator_Apply error : matvec not defined\n");

   if (Op->matvec->ML_id == ML_EXTERNAL)
        Op->matvec->external(Op->data, inlen, din, olen, dout);
   else Op->matvec->internal(Op,       inlen, din, olen, dout);
#ifdef ML_TIMING
   Op->apply_time += (GetClock() - t0);
#endif
   return 0;
}

/* ******************************************************************** */
/* apply the operator to a vector and apply boundary conditions         */
/************************************************************************/

int ML_Operator_ApplyAndResetBdryPts(ML_Operator *Op, int inlen, 
                      double din[], int olen, double dout[])
{
   int i, length, *list;
#ifdef ML_TIMING
   double t0;

   t0 = GetClock();
#endif
   if (Op->matvec->ML_id == ML_EMPTY) 
      pr_error("ML_Operator_ApplyAndRestBdryPts : matvec not defined.\n");

   /* apply grid transfer */

   if (Op->matvec->ML_id == ML_EXTERNAL)
        Op->matvec->external((void*)Op->data, inlen, din, olen, dout);
   else Op->matvec->internal((void*)Op,       inlen, din, olen, dout);

   /* apply boundary condition */

   ML_BdryPts_Get_Dirichlet_Grid_Info(Op->to->BCs, &length, &list);
   for ( i = 0; i < length; i++ ) dout[list[i]] = 0.0;
#ifdef ML_TIMING
   Op->apply_time += (GetClock() - t0);
#endif
   return 0;
}

/* ******************************************************************** */
/* some checking functions                                              */
/*--------------------------------------------------------------------- */

int ML_Operator_Check_Getrow(ML_Operator *Amat, int level, char *str)
{
   int     Nrows, Ncols, i, length, *list;
   double  *t1,*t2,*t3, norm1, norm2;
   ML_Comm *comm;

   if (Amat->getrow->ML_id == ML_EMPTY) return(1);

   comm  = Amat->comm;
   Nrows = Amat->outvec_leng;
   Ncols = Amat->invec_leng;

   if ( Ncols > 0 ) t1 = (double *) ML_allocate(Ncols*sizeof(double) );
   else             t1 = NULL;
   if ( Nrows > 0 ) t2 = (double *) ML_allocate(Nrows*sizeof(double) );
   else             t2 = NULL;
   if ( Nrows > 0 ) t3 = (double *) ML_allocate( Nrows*sizeof(double) );
   else             t3 = NULL;

   for (i = 0; i < Ncols; i++)
      t1[i] = (double) (comm->ML_mypid*2301 + i*i*i*7 + 1);

   if (str[0] == 'R') {
      ML_BdryPts_Get_Dirichlet_Grid_Info(Amat->from->BCs, &length, &list);
      for ( i = 0; i < length; i++ ) t1[list[i]] = 0.0;
      ML_Operator_ApplyAndResetBdryPts(Amat, Ncols, t1, Nrows, t2);
   }
   else ML_Operator_Apply(Amat, Ncols, t1, Nrows, t2);

   norm1 = sqrt(ML_gdot(Nrows, t2, t2, comm));

   ML_getrow_matvec(Amat, t1, Ncols, t3, &Nrows );
   for (i = 0; i < Nrows; i++) t2[i] -= t3[i];
   norm2 = sqrt(ML_gdot(Nrows, t2, t2, comm));
   if (norm2 > norm1*1e-10) {
      norm2 = sqrt(ML_gdot(Nrows, t3, t3, comm));
      if (comm->ML_mypid != 0) return(0);
      printf("Error:\t%s getrow on level %d seems inaccurate\n",str,level);
      printf("\t ||[B] v|| = %e vs. ||B v|| = %e\n",norm2,norm1);
      printf("\twhere [B] v uses %s's getrow routine and B v\n",
	     str);
      printf("\tapplies %s's matrix vector product routine\n",
             str);
   }
   ML_free(t3);
   ML_free(t2);
   ML_free(t1);
   return(0);
}

/* ******************************************************************** */
/* give a name to this operator                                         */
/* ******************************************************************** */

int ML_Operator_Set_Label( ML_Operator *mat, char *label)
{
  int size;

   size = strlen(label) + 1;
   mat->label = (char *) malloc(size*sizeof(char));
   if (mat->label == NULL) pr_error("Not enough space in ML_Operator_Set_Label\n");
   strncpy(mat->label,label,size);
   return(1);
}

/* ******************************************************************** */
/* print this matrix                                                    */
/* ******************************************************************** */

int ML_Operator_Print(ML_Operator *matrix, char label[])
{

   int    i, j;
   int    *bindx;
   double *val;
   int    allocated, row_length;

   if ( matrix->getrow == NULL) return(1);

   allocated = 100;
   bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
   val   = (double *)  ML_allocate( allocated*sizeof(double));

   for (i = 0 ; i < matrix->getrow->Nrows; i++) {
      ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                        &row_length, 0);
      for  (j = 0; j < row_length; j++) {
         printf("%s(%d,%d) = %20.13e;\n",label,i,bindx[j], val[j]);
/*
         printf("%d  %d %20.13e\n",i+1,bindx[j]+1, val[j]);
*/
      }
   }
   ML_free(val);
   ML_free(bindx);
   return 0;
}

/* ******************************************************************** */
/* compute max norm of the matrix                                       */
/* ******************************************************************** */

double ML_Operator_MaxNorm(ML_Operator *matrix, int divide_diag)
{

   int    i, j;
   int    *bindx;
   double *val;
   int    allocated, row_length;
   double sum, largest, diag;

   if ( matrix->getrow == NULL) {
      printf("ML_Operator_MaxNorm: No getrow() function\n");
      return(1.);
   }

   allocated = 100;
   bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
   val   = (double *)  ML_allocate( allocated*sizeof(double));

   largest = 0.;
   for (i = 0 ; i < matrix->getrow->Nrows; i++) {
      ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                        &row_length, 0);
      sum  = 0.;
      diag = 0.;
      for  (j = 0; j < row_length; j++) {
         if (bindx[j] == i) diag = fabs(val[j]);
         sum += fabs(val[j]);
      }
      if (divide_diag == ML_TRUE) {
         if (diag == 0.) printf("ML_Operator_MaxNorm: zero diagonal\n");
         else sum = sum/diag;
      }
      if (sum > largest) largest = sum;
   }
   ML_free(val);
   ML_free(bindx);
   return largest;
}

/* ******************************************************************** */
/* print this matrix (obsolete)                                         */
/* ******************************************************************** */

int ML_Operator_Print_Old(ML_Operator *Amatrix)
{
   int    i, j, m, count, nbytes, Nrows, maxnnz_per_row=500;
   int    *int_buf;
   double *dble_buf;
   int               (*getrowfunc)(void *,int,int*,int,int*,double*,int*);
   ML_GetrowFunc     *getrow_obj;

   /* ------------------------------------------------------------- */
   /* first find out whether the getrow function is available       */
   /* ------------------------------------------------------------- */

   getrow_obj = Amatrix->getrow;
   if ( getrow_obj->ML_id == ML_EXTERNAL) 
   {
      getrowfunc = getrow_obj->external;
   } 
   else 
   {
      getrowfunc = getrow_obj->internal;
   }
   if ( getrowfunc == NULL ) 
   {
      printf("ML_Aggregate_Coarsen Error : no getrow function.\n");
      exit(-1);
   }

   /* ------------------------------------------------------------- */
   /* allocate initial temporary storage space for getrow           */
   /* also allocate space for storing the diagonal (if epsilon>0)   */
   /* ------------------------------------------------------------- */

   Nrows = Amatrix->outvec_leng;
   nbytes = maxnnz_per_row * sizeof( int );
   ML_memory_alloc((void**) &int_buf, nbytes, "OPA");
   nbytes = maxnnz_per_row * sizeof( double );
   ML_memory_alloc((void**) &dble_buf, nbytes, "OPB");

   /* ------------------------------------------------------------- */
   /* find out about how much memory to allocate for the matrix     */
   /* ------------------------------------------------------------- */

   count = 0;
   for ( i = 0; i < Nrows; i++ ) 
   {
      while (getrowfunc(Amatrix->data, 1, &i, maxnnz_per_row, int_buf, 
                          dble_buf, &m) == 0 ) 
      {
         ML_memory_free((void**) &int_buf);
         ML_memory_free((void**) &dble_buf);
         maxnnz_per_row = maxnnz_per_row * 2 + 1; 
         nbytes = maxnnz_per_row * sizeof( int );
         ML_memory_alloc((void**) &int_buf, nbytes, "AGD");
         nbytes = maxnnz_per_row * sizeof( double );
         ML_memory_alloc((void**) &dble_buf,  nbytes, "AGE");
      }
      for ( j = 0; j < m; j++ ) 
         printf("Operator(%d,%d) = %e\n", i+1, int_buf[j]+1, dble_buf[j]);
      count += m;
   }
   ML_memory_free((void**) &int_buf);
   ML_memory_free((void**) &dble_buf);
   return 0;
}

/* ******************************************************************** */
/* Getrow function that is used to drop matrix elements and to collapse */
/* several rows into a block. It is assumed that                        */
/* ML_Operator_AmalgamateAndDropWeak() was previously called to         */
/* properly set up the data structure (data).                           */
/* ******************************************************************** */

int ML_amalg_drop_getrow(void *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   struct amalg_drop *temp;
   int    block_size, row, size, i, j, k, tcol, count;
   int    *tcolumns, tallocated_space;
   double *tvalues, *scaled_diag;
   int offset, status;
   struct ML_GetrowFunc_Struct *amalg_getrow;
   ML_Operator *Amat;
 
   if (N_requested_rows > 1) {
      printf("ML_amalg_drop_getrow: Not implemented for > 1 row at a time\n");
      exit(1);
   }
   temp = (struct amalg_drop *) data;
   Amat = temp->Amat;
   block_size   = temp->block_size;
   amalg_getrow = Amat->getrow;
   scaled_diag  = temp->scaled_diag;

   Amat->data         = temp->original_data;
   Amat->getrow       = temp->original_getrow;
   Amat->invec_leng  *= block_size;
   Amat->outvec_leng *= block_size;

   tallocated_space = allocated_space*block_size*block_size + 1;
   tcolumns     = (int    *) ML_allocate(sizeof(int)*tallocated_space);
   tvalues      = (double *) ML_allocate(sizeof(double)*tallocated_space);
   if (tvalues == NULL) {
      if (tcolumns != NULL) ML_free(tcolumns);
      Amat->data         = temp;
      Amat->getrow       = amalg_getrow;
      Amat->invec_leng  /= block_size;
      Amat->outvec_leng /= block_size;
      return(0);
   }
   offset = 0;
   for (i = 0; i < block_size; i++) {
      row = requested_rows[0]*block_size+i;
      status = ML_Operator_Getrow(Amat, N_requested_rows, &row, 
                                  tallocated_space, &(tcolumns[offset]), 
				  &(tvalues[offset]), &size );
      if (status == 0) {
         ML_free(tvalues); ML_free(tcolumns);
         Amat->data         = temp;
         Amat->getrow       = amalg_getrow;
         Amat->invec_leng  /= block_size;
         Amat->outvec_leng /= block_size;
         return(status);
      }
      if (scaled_diag != NULL) {
         count = 0;
         for (j = offset; j < offset + size; j++) {
            tcol = tcolumns[j];
            if (tvalues[j] != 0.0) {
              if (tvalues[j]*tvalues[j] >= scaled_diag[row]*scaled_diag[tcol]) {
                 tcolumns[offset+count]  = tcolumns[j];
                 tvalues[offset+count++] = tvalues[j];
              }
            }
         }
         size = count;
      }
      tallocated_space -= size;
      offset += size;
   }

   row_lengths[0] = 0;

   for (j = 0; j < offset; j++) {
      tcol = tcolumns[j]/block_size;
      for (k = 0; k < row_lengths[0]; k++) 
         if (tcol == columns[k]) break;

      if (k == row_lengths[0]) {
         if ( allocated_space == row_lengths[0]) {
            ML_free(tvalues); ML_free(tcolumns);
            Amat->data         = temp;
            Amat->getrow       = amalg_getrow;
            Amat->invec_leng  /= block_size;
            Amat->outvec_leng /= block_size;
            return(0);
         }
         values[row_lengths[0]] = 1.;
         columns[row_lengths[0]++] = tcol;
      }
   }
   Amat->data         = temp;
   Amat->getrow       = amalg_getrow;
   Amat->invec_leng  /= block_size;
   Amat->outvec_leng /= block_size;
   ML_free(tvalues); ML_free(tcolumns);
   return(status);
}

/* ******************************************************************** */
/* Restores a matrix that has been modified via                         */
/* ML_Operator_AmalgamateAndDropWeak() back to its original form.       */
/* ******************************************************************** */

int ML_Operator_UnAmalgamateAndDropWeak(ML_Operator *Amat, int block_size,
	double drop_tolerance)
{
   struct amalg_drop *temp;
 
   if ( (block_size > 1) || (drop_tolerance >= 0.0)) {
      temp = (struct amalg_drop *) Amat->data;
      ML_CommInfoOP_Destroy(Amat->getrow->pre_comm);
      ML_memory_free((void**)&(Amat->getrow));
      Amat->data         = temp->original_data;
      Amat->getrow       = temp->original_getrow;
      Amat->invec_leng  *= temp->block_size;
      Amat->outvec_leng *= temp->block_size;
      Amat->num_PDEs     = temp->block_size;
      if (temp->scaled_diag != NULL) free(temp->scaled_diag);
      free(temp);
   }
   return 0;
}
   
/* ******************************************************************** */
/* Modify matrix so that it uses a getrow wrapper that will effectively */
/* drop small values and will collapse several rows into a block row.   */
/* ******************************************************************** */

int ML_Operator_AmalgamateAndDropWeak(ML_Operator *Amat, int block_size, 
               double drop_tolerance)
{
   struct amalg_drop  *new_data;
   int Nneigh, *neighbors, sendleng, rcvleng, *newsend, *newrcv, i, j, k;
   int sendcount, rcvcount, temp, row_length;
   int allocated, *bindx, Nghost, Nrows;
   double *val, *scaled_diag;
   ML_Comm *comm;

   /* create a new widget to hold the amalgamation and drop information */

  comm = Amat->comm;

  if ( (block_size > 1) || (drop_tolerance >= 0.0)) {
     new_data = (struct amalg_drop *) malloc( sizeof(struct amalg_drop) );
     if (new_data == NULL) {
        printf("ML_Operator_AmalgamateAndDropWeak: out of space\n");
        exit(1);
     }
     new_data->original_data   = Amat->data;
     new_data->original_getrow = Amat->getrow;
     new_data->scaled_diag     = NULL;
     new_data->block_size      = block_size;
     new_data->drop_tolerance  = drop_tolerance;
     new_data->Amat            = Amat;
     Nrows                     = Amat->getrow->Nrows;

     /* we need to get the matrix diagonal, scale it by drop_tolerance, */
     /* and store it */

     if ( drop_tolerance >= 0.0) {

        Nneigh    = ML_CommInfoOP_Get_Nneighbors(Amat->getrow->pre_comm);
        neighbors = ML_CommInfoOP_Get_neighbors(Amat->getrow->pre_comm);
        Nghost = 0;
        for (i = 0; i < Nneigh; i++) {
           rcvleng = ML_CommInfoOP_Get_Nrcvlist(Amat->getrow->pre_comm,
                                                neighbors[i]);
           newrcv = ML_CommInfoOP_Get_rcvlist(Amat->getrow->pre_comm, 
                                              neighbors[i]);
           for (j = 0; j < rcvleng; j++) {
              if (newrcv[j] > Nghost + Nrows - 1)
                 Nghost = newrcv[j] - Nrows + 1;
           }
           ML_free(newrcv);
        }
        ML_free(neighbors);

        allocated = 100;
        scaled_diag = (double *) ML_allocate((Nrows+Nghost)*sizeof(double));
        bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
        val   = (double *)  ML_allocate( allocated*sizeof(double));
        if (val == NULL) {
           printf("ML_Operator_AmalgamateAndDropWeak: out of space\n");
           exit(1);
        }

        for (i = 0 ; i < Nrows; i++) {
           ML_get_matrix_row(Amat,1,&i,&allocated,&bindx,&val,&row_length,0);
           for (j = 0; j < row_length; j++) 
              if (bindx[j] == i) break;

           scaled_diag[i] = 0.0;
           if (j != row_length)  scaled_diag[i] = val[j];
           scaled_diag[i] *= drop_tolerance;
           if (scaled_diag[i] < 0.) scaled_diag[i] = -scaled_diag[i];
        }
        ML_free(val);
        ML_free(bindx);
      
        if ( Amat->getrow->pre_comm != NULL )
           ML_exchange_bdry(scaled_diag,Amat->getrow->pre_comm,Nrows, comm, 
                            ML_OVERWRITE);

        new_data->scaled_diag = scaled_diag;
     }

     /* We need to create a new getrow structure */
     /* containing a getrow wrapper              */


     Amat->num_PDEs     = 1;
     Amat->invec_leng  /= block_size;
     Amat->outvec_leng /= block_size;
     Amat->data         = new_data;
     ML_memory_alloc((void**)&(Amat->getrow),sizeof(ML_GetrowFunc),"OF2");
     Amat->getrow->ML_id            = ML_EMPTY;
     Amat->getrow->Nrows            = 0;
     Amat->getrow->pre_comm         = NULL;
     Amat->getrow->post_comm        = NULL;
     Amat->getrow->internal         = NULL;
     Amat->getrow->external         = NULL;
     Amat->getrow->data             = NULL;
     Amat->getrow->use_loc_glob_map = ML_NO;
     Amat->getrow->loc_glob_map     = NULL;
     Amat->getrow->row_map          = NULL;
     ML_Operator_Set_Getrow(Amat, ML_EXTERNAL, 
                            new_data->original_getrow->Nrows/block_size,
                            ML_amalg_drop_getrow);

     /* amalgmation needs a new communication structure. Let's create a new */
     /* communication object and modify it if we are doing amalgmation.     */                   
     ML_CommInfoOP_Clone( &(Amat->getrow->pre_comm),  
                          new_data->original_getrow->pre_comm);

     if (block_size > 1) {
        Nneigh    = ML_CommInfoOP_Get_Nneighbors(Amat->getrow->pre_comm);
        neighbors = ML_CommInfoOP_Get_neighbors(Amat->getrow->pre_comm);

        for (i = 0; i < Nneigh; i++) {
           sendleng = ML_CommInfoOP_Get_Nsendlist(Amat->getrow->pre_comm, 
                                                  neighbors[i]);
           newsend = ML_CommInfoOP_Get_sendlist(Amat->getrow->pre_comm, 
                                                neighbors[i]);
           sendcount = 0;
           for (j = 0 ; j < sendleng; j++) {
              temp = newsend[j]/block_size;

              /* search to see if it is already in the list */
              for (k = 0; k < sendcount; k++) 
                 if ( newsend[k] == temp) break;

              if (k == sendcount) newsend[sendcount++] = temp;
           }
           rcvleng = ML_CommInfoOP_Get_Nrcvlist(Amat->getrow->pre_comm, 
                                                neighbors[i]);
           newrcv = ML_CommInfoOP_Get_rcvlist(Amat->getrow->pre_comm, neighbors[i]);
           rcvcount = 0;
           for (j = 0 ; j < rcvleng; j++) {
              temp = newrcv[j]/block_size;

              /* search to see if it is already in the list */
              for (k = 0; k < rcvcount; k++) 
                 if ( newrcv[k] == temp) break;

              if (k == rcvcount) newrcv[rcvcount++] = temp;
           }
           ML_CommInfoOP_Set_exch_info(Amat->getrow->pre_comm, neighbors[i],
                      rcvcount, newrcv,sendcount, newsend);
           ML_free(newrcv); 
           ML_free(newsend); 
        }
        if (neighbors != NULL) ML_free(neighbors);
     }
  }
  return 0;
}

