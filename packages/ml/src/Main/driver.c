#ifdef ML_MPI
#include <mpi.h>
#include <stdio.h>
#include "ml_include.h"

main(int argc, char *argv[])
{
   int i, j, nrows, *mat_ia, *mat_ja, startRow, mypid, scheme, globalN;
   int nprocs, nullDim;
   double *mat_a, *sol, *rhs;
   MLI_Solver *solver;

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &mypid);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if ( mypid == 0 )
   {
      printf("Test program : no. of processors = %d\n", nprocs);
      printf("Coarsen scheme (1-Uncoupled,2-Coupled) : ");
      scanf("%d", &scheme);
   }
   MPI_Bcast(&scheme, 1, MPI_INT, 0, MPI_COMM_WORLD); 

   solver = MLI_Solver_Create(MPI_COMM_WORLD);
   /*
   nrows = MLI_Solver_Get_IJAFromFile(solver,"matrix.data","rhs.data");
   MLI_Solver_Get_NullSpaceFromFile(solver,"rbm.data");
   */
   globalN = 128 * 128;
   nrows = ML_PDE_GenMat(solver,globalN);
   sol = (double *) malloc( nrows * sizeof(double) );
   for ( i = 0; i < nrows; i++ ) sol[i] = 0.0;
   MLI_Solver_Set_MLNumLevels(solver, 3);
   MLI_Solver_Set_KrylovMethod(solver, 0);
   MLI_Solver_Set_CoarsenScheme(solver,scheme);
   MLI_Solver_Set_StrongThreshold(solver,0.08);
   MLI_Solver_Set_NumPreSmoothings(solver,1);
   MLI_Solver_Set_NumPostSmoothings(solver,1);
   MLI_Solver_Set_PreSmoother(solver,4);
   MLI_Solver_Set_PostSmoother(solver,4);
   MLI_Solver_Set_DampingFactor(solver,1.0);
   MLI_Solver_Setup(solver, sol);
   MLI_Solver_Solve(solver);

   MPI_Finalize();
}
#endif

