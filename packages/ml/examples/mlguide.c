#include <math.h>
#include "ml_include.h"
#ifdef ML_MPI
#include "mpi.h"
#endif
extern int Poisson_getrow(void *A_data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[]);

extern int Poisson_matvec(void *A_data, int in_length, double p[], int out_length,
                   double ap[]);
extern int user_smoothing(void *data, int x_length, double x[],
                   int rhs_length, double rhs[]);


int main(int argc, char *argv[]){


   ML *ml_object;
   int i, N_grids = 3, N_levels;
   double sol[5], rhs[5];
   ML_Aggregate *agg_object;
   ML_Operator *data;

#ifdef ML_MPI
   MPI_Init(&argc,&argv);
#endif
   for (i = 0; i < 5; i++) sol[i] = 0.;
   for (i = 0; i < 5; i++) rhs[i] = 2.;

   ML_Create         (&ml_object, N_grids);

   ML_Init_Amatrix      (ml_object, 0,  5, 5, NULL);
   ML_Set_Amatrix_Getrow(ml_object, 0,  Poisson_getrow, NULL, 5);
   ML_Set_Amatrix_Matvec(ml_object, 0,  Poisson_matvec);

   ML_Aggregate_Create(&agg_object);
   ML_Aggregate_Set_MaxCoarseSize(agg_object,1);
   N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_object, 0,
                                                  ML_INCREASING, agg_object);
   /******** Begin code to set a Jacobi smoother ******

   ML_Gen_Smoother_Jacobi(ml_object, ML_ALL_LEVELS, ML_PRESMOOTHER, 1, ML_DEFAULT);

    ******** End code to set a Jacobi smoother ******/

   /******** Begin code to set a user-defined smoother ******/
   ML_Get_Amatrix(ml_object, 0, &data);
   ML_Set_Smoother(ml_object, 0, ML_PRESMOOTHER, data, user_smoothing,"mine");
   ML_Get_Amatrix(ml_object, 1, &data);
   ML_Set_Smoother(ml_object, 1, ML_PRESMOOTHER, data, user_smoothing,"mine");
   ML_Get_Amatrix(ml_object, 2, &data);
   ML_Set_Smoother(ml_object, 2, ML_PRESMOOTHER, data, user_smoothing,"mine");
   ML_Gen_Solver    (ml_object, ML_MGV, 0, N_levels-1);
   ML_Iterate(ml_object, sol, rhs);
   ML_Aggregate_Destroy(&agg_object);
   ML_Destroy(&ml_object);
   /******** End code to set a user-defined smoother ******/

   printf("answer is %e %e %e %e %e\n",sol[0],sol[1],sol[2],sol[3],sol[4]);

#ifdef ML_MPI
  MPI_Finalize();
#endif
   return(1);
}

int Poisson_getrow(void *A_data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   int count = 0, i, start, row;


   for (i = 0; i < N_requested_rows; i++) {
      if (allocated_space < count+3) return(0);
      start = count;
      row = requested_rows[i];
      if ( (row >= 0) || (row <= (5-1)) ) {
         columns[count] = row; values[count++] = 2.;
         if (row != 0) { columns[count] = row-1; values[count++] = -1.; }
         if (row != (5-1)) { columns[count] = row+1; values[count++] = -1.; }
      }
      row_lengths[i] = count - start;
   }
   return(1);
}
int Poisson_matvec(void *A_data, int in_length, double p[], int out_length,
                   double ap[])
{
   int i;

   for (i = 0; i < 5; i++ ) {
      ap[i] = 2*p[i];
      if (i != 0) ap[i] -= p[i-1];
      if (i != (5-1)) ap[i] -= p[i+1];
   }
   return 0;
}

int user_smoothing(void *data, int x_length, double x[], int rhs_length, double rhs[])
{
   int i;
   double ap[5], omega = .5; /* temp vector and damping factor */
   double *diag;
   ML_Operator *Amat;

   Amat = (ML_Operator *) data;
   ML_Operator_Apply(Amat, x_length, x, rhs_length, ap);
   ML_Operator_Get_Diag(Amat, x_length, &diag);
   
   for (i = 0; i < x_length; i++) x[i] = x[i] + omega*(rhs[i] - ap[i])/diag[i];

   return 0;
}


