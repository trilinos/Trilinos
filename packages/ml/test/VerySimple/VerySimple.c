/* Copied from mlguide.c. 
 * This test user's defined smoothers and getrow's,
 * and the internal Krylov solvers of ML. This is essentially a test
 * on ML itself, without any Trilinos-dependence.
 *
 * MS, 07-Aug-04
 */

#include <math.h>
#include "ml_include.h"
#ifdef ML_MPI
#include "mpi.h"
#endif
extern int Poisson_getrow(ML_Operator *A_data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[]);

extern int Poisson_matvec(ML_Operator *A_data, int in_length, double p[], int out_length,
                   double ap[]);
extern int user_smoothing(ML_Smoother *data, int x_length, double x[],
                   int rhs_length, double rhs[]);


int main(int argc, char *argv[]){


   ML *ml_object;
   int i, N_grids = 3, N_levels;
   double sol[129], rhs[129];
   ML_Aggregate *agg_object;
   ML_Operator *data;
   ML_Krylov *kdata;
   double norm_comp;

#ifdef ML_MPI
   MPI_Init(&argc,&argv);
#endif
   for (i = 0; i < 129; i++) sol[i] = 0.;
   for (i = 0; i < 129; i++) rhs[i] = 2.;
 

   ML_Create         (&ml_object, N_grids);

   ML_Init_Amatrix      (ml_object, 0,  129, 129, NULL);
   ML_Set_Amatrix_Getrow(ml_object, 0,  Poisson_getrow, NULL, 129);
   ML_Set_Amatrix_Matvec(ml_object, 0,  Poisson_matvec);
   ML_Set_PrintLevel(10);

   ML_Aggregate_Create(&agg_object);
   ML_Aggregate_Set_MaxCoarseSize(agg_object,1);
   N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_object, 0,
                                                  ML_INCREASING, agg_object);
   ML_Get_Amatrix(ml_object, 0, &data);
   ML_Set_Smoother(ml_object, 0, ML_BOTH, data, user_smoothing,"mine");
   ML_Get_Amatrix(ml_object, 1, &data);
   ML_Set_Smoother(ml_object, 1, ML_BOTH, data, user_smoothing,"mine");
   ML_Get_Amatrix(ml_object, 2, &data);
   ML_Set_Smoother(ml_object, 2, ML_BOTH, data, user_smoothing,"mine");
   ML_Gen_Solver    (ml_object, ML_MGV, 0, N_levels-1);

   kdata = ML_Krylov_Create(ml_object->comm);
   ML_Krylov_Set_PrintFreq( kdata, 1 );
   ML_Krylov_Set_Method(kdata, ML_CG);
   ML_Krylov_Set_Amatrix(kdata, &(ml_object->Amat[0]));
   ML_Krylov_Set_PreconFunc(kdata, ML_MGVSolve_Wrapper);
   ML_Krylov_Set_Precon(kdata, ml_object);
   ML_Krylov_Set_Tolerance(kdata, 1.e-7);
   ML_Krylov_Solve(kdata, 129, rhs, sol);
   ML_Krylov_Destroy( &kdata );

   ML_Aggregate_Destroy(&agg_object);
   ML_Destroy(&ml_object);

   printf("answer is %e %e %e %e %e\n",sol[0],sol[1],sol[2],sol[3],sol[4]);
   /* subtract what I expect to be the solution */
   sol[0] -= 1.290000e+02;
   sol[1] -= 2.560000e+02;
   sol[2] -= 3.810000e+02;
   sol[3] -= 5.040000e+02;
   sol[4] -= 6.250000e+02;

#ifdef ML_MPI
  MPI_Finalize();
#endif

   norm_comp = 0.0;
   for( i=0 ; i<5 ; ++i ) norm_comp+=sol[i]*sol[i];
   
   if( abs(norm_comp-1.0) > 1e-8 ) return( EXIT_FAILURE );
   else                            return( EXIT_SUCCESS );
}

int Poisson_getrow(ML_Operator *A_data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   int count = 0, i, start, row;



   for (i = 0; i < N_requested_rows; i++) {
      if (allocated_space < count+3) return(0);
      start = count;
      row = requested_rows[i];
      if ( (row >= 0) || (row <= (129-1)) ) {
         columns[count] = row; values[count++] = 2.;
         if (row != 0) { columns[count] = row-1; values[count++] = -1.; }
         if (row != (129-1)) { columns[count] = row+1; values[count++] = -1.; }
      }
      row_lengths[i] = count - start;
   }
   return(1);
}
int Poisson_matvec(ML_Operator *A_data, int in_length, double p[], int out_length,
                   double ap[])
{
   int i;

   for (i = 0; i < 129; i++ ) {
      ap[i] = 2*p[i];
      if (i != 0) ap[i] -= p[i-1];
      if (i != (129-1)) ap[i] -= p[i+1];
   }
   return 0;
}

int user_smoothing(ML_Smoother *data, int x_length, double x[], int rhs_length, double rhs[])
{
   int i;
   double ap[129], omega = .5; /* temp vector and damping factor */
   double *diag;
   ML_Operator *Amat;
   ML_Smoother *smoo;

   smoo    = (ML_Smoother *) data;
   Amat = (ML_Operator *) ML_Get_MySmootherData(smoo);
   ML_Operator_Apply(Amat, x_length, x, rhs_length, ap);
   ML_Operator_Get_Diag(Amat, x_length, &diag);
   
   for (i = 0; i < x_length; i++) x[i] = x[i] + omega*(rhs[i] - ap[i])/diag[i];

   return 0;
}


