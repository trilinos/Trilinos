/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* ******************************************************************** */
/* Functions to create matrix coming out from Poisson problems          */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : February, 2000                                       */
/* ******************************************************************** */
/* ******************************************************************** */

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "ml_struct.h"
#include "ml_pde.h"
#ifdef ML_MPI
#include <mpi.h>

/* ******************************************************************** */
/* generate Poisson matrix in CSR format                                */
/* ******************************************************************** */

int ML_PDE_GenMat(MLI_Solver *solver, int N_nodes) 
{
   int    i, j, k, m, nprocs, nprocs_1d, mypid, mypid_x, mypid_y;
   int    nodeoffset, xoffset, yoffset, index, nbytes;
   int    nnode_1d, nnode_local, nnode_part_1d, **square;
   int    avg_nonzeros_per_row = 5, total_nz, rowcount;
   int    *mat_ia, *mat_ja;
   double alpha=1000.0, diag_sum, *mat_aa, *rhs;
   MPI_Comm comm;

   /* --------------------------------------------------------------- */
   /* box in the processors                                           */
   /* --------------------------------------------------------------- */

   comm   = solver->comm;
   MPI_Comm_rank(comm, &mypid);
   MPI_Comm_size(comm, &nprocs);
   
   nprocs_1d = (int) pow( (double) nprocs, 0.50001 );
   if ( nprocs_1d * nprocs_1d != nprocs ) {
      printf("PDE_GenMat : nprocs should be a square (%d).\n",nprocs_1d);
      exit(1);
   }
   mypid_x = mypid % nprocs_1d;
   mypid_y = mypid / nprocs_1d;

   /* --------------------------------------------------------------- */
   /* process the grid                                                */
   /* --------------------------------------------------------------- */

   nnode_1d = (int) pow( (double) N_nodes, 0.500001 );
   nnode_part_1d = nnode_1d / nprocs_1d;
   if (nnode_part_1d * nprocs_1d != nnode_1d) {
      printf("Error: nnode_part not good. %d %d\n",nnode_part_1d,nnode_1d);
      exit(-1);
   }
   nnode_local = nnode_part_1d * nnode_part_1d;
  
   nbytes = nnode_part_1d * 3 * sizeof(int*);
   ML_memory_alloc((void*) &square, nbytes, "AP1");
   nbytes = nnode_part_1d * 3 * sizeof(int);
   for ( i = 0; i < nnode_part_1d*3; i++ )
      ML_memory_alloc((void*) &(square[i]), nbytes, "AP2");

   for ( j = 0; j < 3; j++ )
   {
      yoffset = j * nnode_part_1d;
      for ( i = 0; i < 3; i++ )
      {
         xoffset = i * nnode_part_1d;
         if ( (mypid_y-1+j) < 0 || (mypid_x-1+i) < 0 ||
              (mypid_y-1+j) >= nprocs_1d || (mypid_x-1+i) >= nprocs_1d )
            nodeoffset = - nnode_local;
         else
            nodeoffset = ((mypid_y-1+j)*nprocs_1d + mypid_x-1+i) * nnode_local;
         for ( k = 0; k < nnode_part_1d; k++ )
            for ( m = 0; m < nnode_part_1d; m++ )
               square[yoffset+k][xoffset+m] = nodeoffset++;            
      }
   }

   total_nz = nnode_local * avg_nonzeros_per_row + 1;
   mat_ia  = (int *) ML_allocate((nnode_local+1) * sizeof(int));
   mat_ja  = (int *) ML_allocate(total_nz * sizeof(int));
   mat_aa  = (double *) ML_allocate(total_nz * sizeof(double));
 
   mat_ia[0] = 0;
   index = 0;
   rowcount = 1;
   xoffset = nnode_part_1d;
   yoffset = nnode_part_1d;
   for ( k = 0; k < nnode_part_1d; k++ )
   {
      for ( m = 0; m < nnode_part_1d; m++ )
      {
         diag_sum = 0.0;
         if ( square[yoffset+k-1][xoffset+m] >= 0 )
         { mat_ja[index] = square[yoffset+k-1][xoffset+m];
           mat_aa[index] = -alpha; index++; diag_sum += alpha;}
         if ( square[yoffset+k][xoffset+m-1] >= 0 )
         { mat_ja[index] = square[yoffset+k][xoffset+m-1];
           mat_aa[index] = -1; index++; diag_sum += 1;}
         mat_ja[index] = square[yoffset+k][xoffset+m];
         mat_aa[index] = 2 * alpha + 2.0;
         index++;
         if ( m < nnode_part_1d-1 || mypid_x < nprocs_1d-1 ) 
         { mat_ja[index] = square[yoffset+k][xoffset+m+1];
           mat_aa[index] = -1; index++; diag_sum += 1;}
         if ( k < nnode_part_1d-1 || mypid_y < nprocs_1d-1 ) 
         { mat_ja[index] = square[yoffset+k+1][xoffset+m];
           mat_aa[index] = -alpha; index++; diag_sum += alpha;}
         mat_ia[rowcount++] = index;
      }
   }

   for ( i = 0; i < nnode_part_1d*3; i++ )
      ML_memory_free((void*) &(square[i]));
   ML_memory_free((void*) &(square));
   rhs  = (double *) ML_allocate(nnode_local * sizeof(double));
   for ( i = 0; i < nnode_local; i++ ) rhs[i] = 1.0; 

   solver->nRows = nnode_local;
   solver->mat_ia = mat_ia;
   solver->mat_ja = mat_ja;
   solver->mat_a  = mat_aa;
   solver->rhs    = rhs;
   return nnode_local;
}

#endif
