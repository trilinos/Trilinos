/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for the ML_Krylov structure                                */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : December, 1999                                       */
/* ******************************************************************** */

#include <stdlib.h>
#include "ml_krylov.h"

/* ******************************************************************** */
/* Create an ML_Krylov and initialize relevant fields.                  */
/* ******************************************************************** */

ML_Krylov *ML_Krylov_Create(ML_Comm *comm)
{
   ML_Krylov *temp;

   ML_memory_alloc((void**) &temp, sizeof(ML_Krylov), "KR1");
   temp->ML_id             = ML_ID_KRYLOVDATA;
   temp->ML_method         = 1;
   temp->ML_gmres_dim      = 300;
   temp->ML_cgstabl_dim    = 2;
   temp->ML_max_iterations = 1000;
   temp->ML_tolerance      = 1.0e-6;
   temp->ML_matrix         = NULL;
   temp->ML_com            = comm;
   temp->ML_matrix         = NULL;
   temp->ML_precon         = NULL;
   temp->ML_precfcn        = NULL;
   temp->ML_eigen          = 0;
   temp->ML_nonsym_eigen   = 0;
   temp->ML_eigen_max      = 0.0;
   temp->diag_scale        = NULL;
   temp->ML_print_freq     = 1;
   temp->ML_dont_scale_by_diag = 0; /* don't scale by diagonal */
                                    /* when computing eigenvalues */

   return(temp);
}

/* ******************************************************************** */
/* destructor                                                           */
/* ******************************************************************** */

int ML_Krylov_Destroy( ML_Krylov **data)
{
   if ( (*data)->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Destroy error : wrong object.\n");
      exit(-1);
   }
   if ( (*data)->diag_scale != NULL ) ML_free((*data)->diag_scale);
   ML_memory_free((void**) data);
   return 0;
}

/* ******************************************************************** */
/* get the communicator                                                 */
/* ******************************************************************** */

ML_Comm * ML_Krylov_Get_Comm(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Get_Comm error : wrong object.\n");
      exit(-1);
   }
   return data->ML_com;
}

/* ******************************************************************** */
/* choose a Krylov method                                               */
/* ******************************************************************** */

int ML_Krylov_Set_Method(ML_Krylov *data, int method)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_Method error : wrong object.\n");
      exit(-1);
   }
   if      ( method == 0 ) data->ML_method = 0;
   else if ( method == 1 ) data->ML_method = 1;
   else                    data->ML_method = 1;
   return 0;
}

/* ******************************************************************** */
/* set print frequency                                                  */
/* ******************************************************************** */

int ML_Krylov_Set_PrintFreq(ML_Krylov *data, int n)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_PrintFreq error : wrong object.\n");
      exit(-1);
   }
   data->ML_print_freq = n;
   return 0;
}

/* ******************************************************************** */
/* get print frequency                                                  */
/* ******************************************************************** */

int ML_Krylov_Get_PrintFreq(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Get_PrintFreq error : wrong object.\n");
      exit(-1);
   }
   return data->ML_print_freq;
}

/* ******************************************************************** */
/* set the dimension of Krylov subspace to keep in restarted GMRES      */
/* ******************************************************************** */

int ML_Krylov_Set_GMRESSize(ML_Krylov *data, int size)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_GMRESSize error : wrong object.\n");
      exit(-1);
   }
   if ( size > 1 ) data->ML_gmres_dim = size;
   else            data->ML_gmres_dim = 20;
   return 0;
}

/* ******************************************************************** */
/* get the dimension of Krylov subspace to keep in restarted GMRES      */
/* ******************************************************************** */

int ML_Krylov_Get_GMRESSize(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Get_GMRESSize error : wrong object.\n");
      exit(-1);
   }
   return data->ML_gmres_dim;
}

/* ******************************************************************** */
/* set the dimension of Krylov subspace to keep in BiCGSTABL            */
/* ******************************************************************** */

int ML_Krylov_Set_BICGSTABLSize(ML_Krylov *data, int size)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_BICGSTABLSize error : wrong object.\n");
      exit(-1);
   }
   if      ( size == 1 ) data->ML_cgstabl_dim = size;
   else if ( size == 2 ) data->ML_cgstabl_dim = size;
   else if ( size == 4 ) data->ML_cgstabl_dim = size;
   else                  data->ML_cgstabl_dim = 2;
   return 0;
}

/* ******************************************************************** */
/* get the dimension of Krylov subspace to keep in BICGSTABL            */
/* ******************************************************************** */

int ML_Krylov_Get_BICGSTABLSize(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Get_BICGSTABLSize error : wrong object.\n");
      exit(-1);
   }
   return data->ML_cgstabl_dim;
}

/* ******************************************************************** */
/* set the matrix                                                       */
/* ******************************************************************** */

int ML_Krylov_Set_Amatrix(ML_Krylov *data, ML_Operator *mat)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_Amatrix error : wrong object.\n");
      exit(-1);
   }
   data->ML_matrix = mat;
   return 0;
}

/* ******************************************************************** */
/* set the matrix                                                       */
/* ******************************************************************** */

ML_Operator* ML_Krylov_Get_Amatrix(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Get_Amatrix error : wrong object.\n");
      exit(-1);
   }
   return data->ML_matrix;
}

/* ******************************************************************** */
/* set maximum number of iterations                                     */
/* ******************************************************************** */

int ML_Krylov_Set_MaxIterations(ML_Krylov *data, int iter )
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_MaxIterations error : wrong object.\n");
      exit(-1);
   }
   data->ML_max_iterations = iter;
   return 0;
}

/* ******************************************************************** */
/* get maximum number of iterations                                     */
/* ******************************************************************** */

int ML_Krylov_Get_MaxIterations(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Get_MaxIterations error : wrong object.\n");
      exit(-1);
   }
   return data->ML_max_iterations;
}

/* ******************************************************************** */
/* set tolerance for convergence                                        */
/* ******************************************************************** */

int ML_Krylov_Set_Tolerance(ML_Krylov *data, double tol )
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_Tolerance error : wrong object.\n");
      exit(-1);
   }
   data->ML_tolerance = tol;
   return 0;
}

/* ******************************************************************** */
/* set diagonal inverse                                                 */
/* ******************************************************************** */

int ML_Krylov_Set_Diagonal(ML_Krylov *data, int leng, double *diag )
{
   int  i;

   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_Diagonal error : wrong object.\n");
      exit(-1);
   }
   if ( leng > 0 ) data->diag_scale = (double*) ML_allocate(leng * sizeof(double));
   else            data->diag_scale = NULL;
   printf("set diag = %d\n", leng);

   for ( i = 0; i < leng; i++ ) data->diag_scale[i] = diag[i];
   return 0;
}

/* ******************************************************************** */
/* get tolerance for convergence                                        */
/* ******************************************************************** */

double ML_Krylov_Get_Tolerance(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Get_Tolerance error : wrong object.\n");
      exit(-1);
   }
   return data->ML_tolerance;
}

/* ******************************************************************** */
/* set preconditioner object                                            */
/* ******************************************************************** */

int ML_Krylov_Set_Precon(ML_Krylov *data, void *prec)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_Precon error : wrong object.\n");
      exit(-1);
   }
   data->ML_precon = prec;
   return 1;
}

/* ******************************************************************** */
/* get preconditioner object                                            */
/* ******************************************************************** */

void *ML_Krylov_Get_Precon(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Get_Precon error : wrong object.\n");
      exit(-1);
   }
   return data->ML_precon;
}

/* ******************************************************************** */
/* set preconditioner function                                          */
/* ******************************************************************** */

int ML_Krylov_Set_PreconFunc(ML_Krylov *data, 
                            int (*func)(void*,int,double*,int,double*))
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_PreconFunc error : wrong object.\n");
      exit(-1);
   }
   data->ML_precfcn = func;
   return 1;
}

/* ******************************************************************** */
/* solve the linear system                                              */
/* ******************************************************************** */
#include "ml_cg.h"
#include "ml_gmres.h"
#include "ml_bicgstabl.h"

int ML_Krylov_Solve(ML_Krylov *data,int leng,double *invec,double* outvec)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Solve error : wrong object.\n");
      exit(-1);
   }

   if ( data->ML_eigen == 1 )
   {
#ifdef RST_MODIF
      if (data->ML_dont_scale_by_diag == 1)
	ML_CG_ComputeEigenvalues(data, leng, ML_FALSE); 
      else ML_CG_ComputeEigenvalues(data, leng, ML_TRUE); 
#else
#ifndef MB_MODIF
      if (data->ML_dont_scale_by_diag == 1)
	ML_CG_ComputeEigenvalues(data, leng, ML_FALSE); 
      else ML_CG_ComputeEigenvalues(data, leng, ML_TRUE); 

#else
      /* no diag. in smoother */
      ML_CG_ComputeEigenvalues(data, leng, ML_FALSE);
#endif
#endif
      data->ML_eigen = 0;
      return 0;
   }

   else if ( data->ML_nonsym_eigen == 1 )
   {
#ifdef RST_MODIF
      if (data->ML_dont_scale_by_diag == 1)
	ML_Power_ComputeEigenvalues(data, leng, ML_FALSE); 
      else ML_Power_ComputeEigenvalues(data, leng, ML_TRUE); 
#else
#ifndef MB_MODIF
      if (data->ML_dont_scale_by_diag == 1)
	ML_Power_ComputeEigenvalues(data, leng, ML_FALSE); 
      else ML_Power_ComputeEigenvalues(data, leng, ML_TRUE); 
#else
      /* no diag. in smoother */
      ML_Power_ComputeEigenvalues(data, leng, ML_FALSE); 
#endif
#endif
      data->ML_nonsym_eigen = 0;
      return 0;
   }

   switch ( data->ML_method )
   {
      case ML_CG :   ML_CG_Solve(data, leng, invec, outvec);
                break;
      case ML_GMRES : ML_GMRES_Solve(data, leng, invec, outvec);
                break;
      default : ML_BICGSTABL_Solve(data, leng, invec, outvec);
                break;
   }
   return 0;
}

/* ******************************************************************** */
/* wrapper for diagonal scale                                           */
/* ******************************************************************** */

int ML_DiagScale_Wrapper(void *data,int leng,double *outvec,int leng2,
                         double* invec)
{
   int       i;
   double    *diag_scale;
   ML_Krylov *kry_obj;

   kry_obj = (ML_Krylov *) data;

   if ( leng != leng2 )
   {
      printf("ML_DiagScale_Wrapper ERROR : lengths do not match.\n");
      exit(0);
   }
   diag_scale = kry_obj->diag_scale;

   for ( i = 0; i < leng; i++ ) outvec[i] = invec[i] * diag_scale[i];

   return 0;
}

/* ******************************************************************** */
/* wrapper for AMG solve                                                */
/* ******************************************************************** */

int ML_MGVSolve_Wrapper(void *data,int leng,double *outvec,int leng2,
                         double* invec)
{
   ML *ml_ptr;

   ml_ptr = (ML *) data;

   if ( leng != leng2 )
   {
      printf("ML_DiagScale_Wrapper ERROR : lengths do not match.\n");
      exit(0);
   }
   ML_Solve_MGV(ml_ptr, invec, outvec);

   return 0;
}

/* ******************************************************************** */
/* wrapper for AMG solve                                                */
/* ******************************************************************** */

int ML_AMGVSolve_Wrapper(void *data,int leng,double *outvec,int leng2,
                         double* invec)
{
   ML *ml_ptr;

   ml_ptr = (ML *) data;

   if ( leng != leng2 )
   {
      printf("ML_DiagScale_Wrapper ERROR : lengths do not match.\n");
      exit(0);
   }
   ML_Solve_AMGV(ml_ptr, invec, outvec);

   return 0;
}

/* ******************************************************************** */
/* set to compute max eigenvalue                                        */
/* ******************************************************************** */

int ML_Krylov_Set_ComputeEigenvalues(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_ComputeEigenvalues error : wrong object.\n");
      exit(-1);
   }
   data->ML_eigen = 1;
   return 0;
}
int ML_Krylov_Set_ComputeNonSymEigenvalues(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Set_ComputeNonSymEigenvalues error : wrong object.\n");
      exit(-1);
   }
   data->ML_nonsym_eigen = 1;
   return 0;
}

/* ******************************************************************** */
/* Indicate whether the matrix should be scaled by the diagonal         */
/* when computing eigenvalues.                                          */
/* ******************************************************************** */
int ML_Krylov_Set_DiagScaling_Eig(ML_Krylov *data, int scale)
{
  /* scale = 1: use diagonal scaling       */
  /* scale = 0: don't use diagonal scaling */

  if (scale == 1)
    data->ML_dont_scale_by_diag = 0;
  else if (scale == 0) 
    data->ML_dont_scale_by_diag = 1;
  else {
    printf("ML_Krylov_Set_DiagScaling_Eig: Unknown scaling option %d\n",
	   scale);
    return 1;
  }
  return 0;
}

/* ******************************************************************** */
/* get to spectral radius                                               */
/* ******************************************************************** */

double ML_Krylov_Get_MaxEigenvalue(ML_Krylov *data)
{
   if ( data->ML_id != ML_ID_KRYLOVDATA ) 
   {
      printf("ML_Krylov_Get_MaxEigenvalue error : wrong object.\n");
      exit(-1);
   }
   return data->ML_eigen_max;
}

