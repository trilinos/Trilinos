#include <math.h>
#include "ml_include.h"
#ifdef ML_MPI
#include "mpi.h"
#endif

extern int Poisson_getrow(ML_Operator *mat_in, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[]);

extern int Poisson_matvec(ML_Operator *mat_in, int in_length, double p[], int out_length,
                   double ap[]);
extern int Poisson_comm(double x[], void *A_data);
extern int send_msg(char *send_buffer,  int length, int neighbor);
extern int recv_msg(char *recv_buffer,  int length, int neighbor, USR_REQ *request);
extern int post_msg(char *recv_buffer,  int length, int neighbor, USR_REQ *request);


/* Simple example corresponding to Poisson on a line with a total of 5 grid points */
/* Thus, the global matrix is                                                      */
/*                                                                                 */
/*                           |  2  -1              |                               */
/*                           | -1   2  -1          |                               */
/*                           |     -1   2  -1      |                               */
/*                           |         -1   2  -1  |                               */
/*                           |             -1   2  |                               */
/*                                                                                 */
/* Processor 0 is assigned global rows 0 and 4 which are stored locally as 1 and   */
/* 0. Processor 1 is assigned global rows 1-3 which are stored locally as 0-2      */
/* Processor 0's local matrix is                                                   */
/*                                                                                 */
/*                           |  2          -1      |                               */
/*                           |      2   1      -1  |                               */
/*                                                                                 */
/* and processor 1's local matrix is                                               */
/*                                                                                 */
/*                           |  2  -1          -1  |                               */
/*                           | -1   2  -1          |                               */
/*                           |     -1   2  -1      |                               */
/*                                                                                 */
/* See the ML guide for more details.                                              */

int main(int argc, char *argv[]){


   ML *ml_object;
   int i, N_grids = 3, N_levels;
   double sol[5], rhs[5];
   ML_Aggregate *agg_object;
   int proc, nlocal, nlocal_allcolumns;

   MPI_Init(&argc,&argv);
   ML_Set_PrintLevel(15);

   for (i = 0; i < 5; i++) sol[i] = 0.;
   for (i = 0; i < 5; i++) rhs[i] = 2.;


   ML_Create         (&ml_object, N_grids);
   proc = ml_object->comm->ML_mypid;
   if (ml_object->comm->ML_nprocs != 2) {
      if (proc == 0) printf("Must be run on two processors\n");
      ML_Destroy(&ml_object);
      MPI_Finalize();
      exit(1);
   }

   if     (proc == 0) {nlocal = 2; nlocal_allcolumns = 4;}
   else if (proc == 1){nlocal = 3; nlocal_allcolumns = 5;}
   else               {nlocal = 0; nlocal_allcolumns = 0;}

   ML_Init_Amatrix      (ml_object, 0,  nlocal, nlocal, &proc);
   ML_Set_Amatrix_Getrow(ml_object, 0,  Poisson_getrow, Poisson_comm,
                         nlocal_allcolumns);
   ML_Set_Amatrix_Matvec(ml_object, 0,  Poisson_matvec);

   ML_Aggregate_Create(&agg_object);
   ML_Aggregate_Set_MaxCoarseSize(agg_object,1);
   N_levels = ML_Gen_MGHierarchy_UsingAggregation(ml_object, 0,
                                                  ML_INCREASING, agg_object);
   ML_Gen_Smoother_Jacobi(ml_object, ML_ALL_LEVELS, ML_PRESMOOTHER, 1, ML_DEFAULT);

   ML_Gen_Solver    (ml_object, ML_MGV, 0, N_levels-1);
   ML_Iterate(ml_object, sol, rhs);

   if (proc == 0) {
      printf("sol(0) = %e\n",sol[1]);
      fflush(stdout);
   }
   ML_Comm_GsumInt(ml_object->comm,1);    /* just used for synchronization */
   if (proc == 1) {
      printf("sol(1) = %e\n",sol[0]);
      printf("sol(2) = %e\n",sol[1]);
      printf("sol(3) = %e\n",sol[2]);
      fflush(stdout);
   }
   ML_Comm_GsumInt(ml_object->comm,1);    /* just used for synchronization */
   if (proc == 0) {
      printf("sol(4) = %e\n",sol[0]);
      fflush(stdout);
   }
   ML_Aggregate_Destroy(&agg_object);
   ML_Destroy(&ml_object);

   MPI_Finalize();

   return 0;
}

/* Application specific getrow. */
int Poisson_getrow(ML_Operator *mat_in, int N_requested_rows, int requested_rows[],
   int allocated_space, int cols[], double values[], int row_lengths[])
{
   int m = 0, i, row, proc, *itemp, start;

   itemp  = (int *) ML_Get_MyGetrowData(mat_in);

   proc  = *itemp;

   for (i = 0; i < N_requested_rows; i++) {
      row = requested_rows[i];
      if (allocated_space < m+3) return(0);
      values[m] = 2; values[m+1] = -1; values[m+2] = -1;
      start = m;
      if (proc == 0) {
         if (row == 0) {cols[m++] = 0; cols[m++] = 2;               }
         if (row == 1) {cols[m++] = 1; cols[m++] = 3;}
      }
      if (proc == 1) {
         if (row == 0) {cols[m++] = 0; cols[m++] = 1; cols[m++] = 4;}
         if (row == 1) {cols[m++] = 1; cols[m++] = 0; cols[m++] = 2;}
         if (row == 2) {cols[m++] = 2; cols[m++] = 1; cols[m++] = 3;}
      }
      row_lengths[i] = m - start;
   }
   return(1);
}

/* Application specific matrix-vector product. */
int Poisson_matvec(ML_Operator *mat_in, int in_length, double p[], int out_length,
                   double ap[])
{
   int i, proc, *itemp;
   double new_p[5];

   itemp = (int *) ML_Get_MyMatvecData(mat_in);
   
   proc  = *itemp; 

   for (i = 0; i < in_length; i++) new_p[i] = p[i];
   Poisson_comm(new_p, &proc);

   for (i = 0; i < out_length; i++) ap[i] = 2.*new_p[i];

   if (proc == 0) {
      ap[0] -= new_p[2];
      ap[1] -= new_p[3];
   }
   if (proc == 1) {
      ap[0] -= new_p[1]; ap[0] -= new_p[4];
      ap[1] -= new_p[2]; ap[1] -= new_p[0];
      ap[2] -= new_p[3]; ap[2] -= new_p[1];
   }
   return 0;
}

/* Communication routine that should be performed before doing a matvec */
/* See ML Guide.                                                        */
int Poisson_comm(double x[], void *A_data)
{
   int    proc, neighbor, length, *itemp;
   double send_buffer[2], recv_buffer[2];
   MPI_Request request;

   itemp = (int *) A_data;
   proc  = *itemp; 

   length = 2;
   if (proc == 0) {
      neighbor = 1;
      send_buffer[0] = x[0]; send_buffer[1] = x[1];
      post_msg((char *) recv_buffer,  length, neighbor, &request);
      send_msg((char *) send_buffer,  length, neighbor);
      recv_msg((char *) recv_buffer,  length, neighbor, &request);
      x[2] = recv_buffer[1]; x[3] = recv_buffer[0];
   }
   else {
      neighbor = 0;
      send_buffer[0] = x[0]; send_buffer[1] = x[2];
      post_msg((char *) recv_buffer,  length, neighbor, &request);
      send_msg((char *) send_buffer,  length, neighbor);
      recv_msg((char *) recv_buffer,  length, neighbor, &request);
      x[3] = recv_buffer[1]; x[4] = recv_buffer[0];
   }
   return 0;
}

/* Simple communication wrappers for use with MPI */

int send_msg(char *send_buffer,  int length, int neighbor)
{
   ML_Comm_Send(send_buffer, length*sizeof(double), neighbor, 123,
                 MPI_COMM_WORLD);
   return 0;
}
int recv_msg(char *recv_buffer,  int length, int neighbor, USR_REQ *request)
{
   MPI_Status status;
   MPI_Wait(request, &status);
   return 0;
}
int post_msg(char *recv_buffer,  int length, int neighbor, USR_REQ *request)
{
   int type = 123;
   ML_Comm_Irecv(recv_buffer, length*sizeof(double), &neighbor,
                  &type, MPI_COMM_WORLD, request);
   return 0;
}
