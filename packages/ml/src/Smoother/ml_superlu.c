/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* -- SuperLU routine (version 1.1) --                                  */
/* Univ. of California Berkeley, Xerox Palo Alto Research Center,       */
/* and Lawrence Berkeley National Lab.                                  */
/* November 15, 1997                                                    */
/* ******************************************************************** */
/* ML interface to SuperLU :                                            */
/*   Author : Charles Tong (LLNL) and Ray Tuminaro (SNL)                */
/*   Date   : April, 1998                                               */
/* ******************************************************************** */
#include <stdlib.h>
#include "ml_common.h"
#ifdef SUPERLU
#include "dsp_defs.h"
#include "util.h"
#else
#ifdef DSUPERLU
#include "mpi.h"
#include <malloc.h>
#include "superlu_ddefs.h"
#endif
#endif

#include "ml_comm.h"
#include "ml_struct.h"
#include "ml_solver.h"
#include "ml_superlu.h"
#define REPLACE 1
#define NEQU    0

/* This macro should be defined to reconcile name changes between
   SuperLU 1.0 & SuperLU 2.0.*/
#ifdef ML_SUPERLU2_0
#define DN SLU_DN
#define _D SLU_D
#define GE SLU_GE
#define NR SLU_NR
#endif

/* ******************************************************************** */
/* This subroutine calls the SuperLU subroutine to perform LU           */
/* factorization of a given matrix                                      */
/* ******************************************************************** */

int ML_SuperLU_Solve(ML_Solver *vsolver,int ilen,double *x,int olen,double *rhs)
{
#ifdef SUPERLU
   int            i, n, info, flag, *perm_r, *perm_c, permc_spec;
   int            N_local, offset, *etree, panel_size, lwork;
   double         *local_x, *local_rhs, *R, *C, *ferr, *berr;
   double         rpg, rcond;
   char           fact[1], equed[1], trans[1], refact[1];
   void           *work=NULL;
   ML_Comm        *comm;
   ML_Solver      *solver;
   SuperMatrix    *A, *L, *U, B, X;
   factor_param_t iparam;
   mem_usage_t    mem_usage;

#ifdef ML_TIMING_DETAILED
   double         t0;
#endif

   /* ------------------------------------------------------------- */
   /* fetch the sparse matrix and other parameters                  */
   /* ------------------------------------------------------------- */

   if ( ilen != olen ) 
   {
      printf("ML_SuperLU_Solve error : lengths not matched.\n");
      exit(1);
   }
   solver   = (ML_Solver *) vsolver;

   comm     = (ML_Comm *) solver->void_params1;
   N_local  = (int) solver->dble_params1[0];
   offset   = (int) solver->dble_params1[1];
   n        = (int) solver->dble_params1[2];
   flag     = solver->reuse_flag;

   /* ------------------------------------------------------------- */
   /* if factorization has not been done, allocate space for it     */
   /* ------------------------------------------------------------- */

   if ( flag == -999 ) 
   {
      A = (SuperMatrix *) solver->Mat1;
      if ( A != NULL )
      {
         ML_memory_free((void**)&(((NRformat *) A->Store)->colind));
         ML_memory_free((void**)&(((NRformat *) A->Store)->rowptr));
         ML_memory_free((void**)&(((NRformat *) A->Store)->nzval));
         SUPERLU_FREE( ((NRformat *) A->Store)->nzval);
         SUPERLU_FREE( A->Store );
         ML_memory_free((void**) &A);
      }
      solver->Mat1 = NULL;
      L = (SuperMatrix *) solver->Mat2;
      if ( L != NULL ) 
      {
         Destroy_SuperNode_Matrix(L);
         ML_memory_free((void**) &L);
      }
      solver->Mat2 = NULL;
      U = (SuperMatrix *) solver->Mat3;
      if ( U != NULL )
      {
         SUPERLU_FREE( ((NRformat *) U->Store)->colind);
         SUPERLU_FREE( ((NRformat *) U->Store)->rowptr);
         SUPERLU_FREE( ((NRformat *) U->Store)->nzval);
         SUPERLU_FREE( U->Store );
         ML_memory_free((void**) &U);
      }
      solver->Mat3 = NULL;
      perm_r = (int *) solver->int_params1;
      if (perm_r != NULL) ML_memory_free((void**) &(solver->int_params1));
      solver->int_params1 = NULL;
      perm_c = (int *) solver->int_params2;
      if (perm_c != NULL) ML_memory_free((void**) &(solver->int_params2));
      solver->int_params2 = NULL;
#ifdef ML_TIMING_DETAILED
      if ( comm->ML_mypid == 0 )
         printf("Total SuperLU solve time = %e\n", solver->dble_data);
#endif
      return 0;
   } 
   else if ( flag == 0 ) 
   {
      A = (SuperMatrix *) solver->Mat1;
      ML_memory_alloc((void**) &perm_c, n * sizeof(int), "LU2" );
      ML_memory_alloc((void**) &perm_r, 2 * n * sizeof(int), "LU3" );
      solver->int_params1 = perm_r;
      solver->int_params2 = perm_c;
      permc_spec = 0;
      get_perm_c(permc_spec, A, perm_c);
      ML_memory_alloc((void**) &L, sizeof(SuperMatrix), "LU4" );
      ML_memory_alloc((void**) &U, sizeof(SuperMatrix), "LU5" );
      solver->Mat2 = (void *) L;
      solver->Mat3 = (void *) U;
   } 
   else 
   {
      A = (SuperMatrix *) solver->Mat1;
      L = (SuperMatrix *) solver->Mat2;
      U = (SuperMatrix *) solver->Mat3;
      perm_r = (int *) solver->int_params1;
      perm_c = (int *) solver->int_params2;
   }

   /* ------------------------------------------------------------- */
   /* gather from all processors the complete right hand side       */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &local_rhs, n*sizeof(double),"LU1" );
   ML_memory_alloc((void**) &local_x,   n*sizeof(double),"LU2" );
   for ( i = 0; i < N_local; i++ ) local_rhs[i] = rhs[i];
   i = N_local;
   ML_Comm_GappendDouble((ML_Comm *) comm, local_rhs, &i, n);

   /* ------------------------------------------------------------- */
   /* create right hand side vector conforming to SuperLU format    */
   /* ------------------------------------------------------------- */

   dCreate_Dense_Matrix(&B, n, 1, local_rhs, n, DN, _D, GE);
   dCreate_Dense_Matrix(&X, n, 1, local_x,   n, DN, _D, GE);

   /* ------------------------------------------------------------- */
   /* perform LU decomposition and then solve                       */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &etree,   n*sizeof(int),"LU3" );
   panel_size               = sp_ienv(1);
   iparam.panel_size        = panel_size;
   iparam.relax             = sp_ienv(2);
   iparam.diag_pivot_thresh = 1.0;
   iparam.drop_tol          = -1;
   lwork                    = 0;
   if ( flag == 0 )   *fact = 'N';
   else               *fact = 'F';
   *equed                   = 'N';
   *trans                   = 'N';
   *refact                  = 'N';
   R    = (double *) SUPERLU_MALLOC(n * sizeof(double));
   C    = (double *) SUPERLU_MALLOC(n * sizeof(double));
   ferr = (double *) SUPERLU_MALLOC(sizeof(double));
   berr = (double *) SUPERLU_MALLOC(sizeof(double));

#ifdef ML_TIMING_DETAILED
   t0 = GetClock();
#endif

   dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
          equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
          ferr, berr, &mem_usage, &info);

#ifdef ML_TIMING_DETAILED
   t0 = GetClock() - t0;
   solver->dble_data += t0;
#endif

   if ( info != 0 && info != n+1 )
   {
      printf("ML_SuperLU_Solve : error coming from dgssvx %d\n", info);
      exit(1);
   } 
   else if ( solver->reuse_flag == 0 )
   {
#ifdef ML_DEBUG_SUPERLU
      if ( rcond != 0.0 && offset == 0 )
         printf("ML_SuperLU_Solve : condition number = %e\n", 1.0/rcond);
      else if ( offset == 0 )
         printf("ML_SuperLU_Solve : Recip. condition number = %e\n", rcond);
#endif
   }

   /* ------------------------------------------------------------- */
   /* extract the local solution sub-vector and then clean up       */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < N_local; i++ ) x[i] = local_x[i+offset];

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   SUPERLU_FREE (R);
   SUPERLU_FREE (C);
   SUPERLU_FREE (ferr);
   SUPERLU_FREE (berr);
   ML_memory_free( (void **) &local_x );
   ML_memory_free( (void **) &local_rhs );
   ML_memory_free( (void **) &etree );
   Destroy_SuperMatrix_Store(&B);
   Destroy_SuperMatrix_Store(&X);
   solver->reuse_flag = 1;
#else
#ifdef DSUPERLU
   int                flag, N_local, offset;
   double            *local_rhs;
   ML_Comm           *mlcomm;
   ML_Solver         *solver;
   superlu_options_t  options;
   superlu_options_t *optionsptr;
   SuperMatrix       *A;
   ScalePermstruct_t *ScalePermstruct;
   SuperLUStat_t      stat;
   LUstruct_t        *LUstruct;
   gridinfo_t        *mygrid;
   ML_Lugrid         *lugrid_tiles;
   double             berr;
   int_t              i, n, nprow, npcol, nprocs; /* short integers */
   int                iam, info, nrhs, color, key;
   int                q, g, l, k, mygroup;
   int                stile, mtile, ltile, ntile, tsz, tsz2;
   int_t             *usermap;
   /* In-lining superlu_gridmap */
   MPI_Group mpi_base_group, superlu_grp;
   int                mycol, myrow, j;
   /* heap info arguments */
   int fragments, total_free, largest_free, total_used;

   /* ------------------------------------------------------------- */
   /* fetch the sparse matrix and other parameters                  */
   /* ------------------------------------------------------------- */

   if ( ilen != olen ) 
   {
      printf("ML_SuperLU_Solve error : lengths not matched.\n");
      exit(1);
   }
   solver          = (ML_Solver *) vsolver;
   A               = (SuperMatrix *) solver->Mat1;
   mlcomm          = (ML_Comm *) solver->void_params1;
   iam             = mlcomm->ML_mypid;
   nprocs          = mlcomm->ML_nprocs;
   N_local         = (int) solver->dble_params1[0];
   offset          = (int) solver->dble_params1[1];
   i               = (int) solver->dble_params1[2];
   n               = i;   /* n is a short integer */
   flag            = solver->reuse_flag;
   optionsptr      = &options;
   LUstruct        = (LUstruct_t *) solver->LUspl;
   ScalePermstruct = (ScalePermstruct_t *) solver->PERMspl;

   /* ------------------------------------------------------------- */
   /* if factorization has not been done, allocate space for it     */
   /* Fetching the factorization (flag=1) is unneccessary           */
   /* ------------------------------------------------------------- */

   if ( flag == -999 ) 
   {
     if( iam == 0 )printf("ml_superlu: clean up\n"); /* dmd */

     /* deallocate storage and clean up */
     info = flag;
     if ( A != NULL ) 
     {
        ML_memory_free((void*)&(((NCformat *) A->Store)->rowind));
        ML_memory_free((void*)&(((NCformat *) A->Store)->colptr));
        ML_memory_free((void*)&(((NCformat *) A->Store)->nzval));
        SUPERLU_FREE( A->Store );
        ML_memory_free((void**) &A);
     }
     solver->Mat1 = NULL;
     mygroup = (int ) solver->ML_subgroup;
     lugrid_tiles    = solver->gridtiles;
     mygrid = &((lugrid_tiles[mygroup]).grid);
     Destroy_LU(n, mygrid, LUstruct);
     ScalePermstructFree(ScalePermstruct);
     LUstructFree(LUstruct);
     ML_free(LUstruct);
     ML_free(ScalePermstruct);
     solver->PERMspl = NULL;
     solver->LUspl = NULL;
     superlu_gridexit(mygrid);
     ML_free(lugrid_tiles);
     solver->gridtiles = NULL;
     return 0;
   } 
   else if ( flag == 0 ) 
   {
      ML_SuperLU_Set_Tile(nprocs, &tsz ,&stile, &mtile, &ltile);
      ntile = stile + mtile + ltile;
      tsz2 = tsz * tsz;
      usermap = (int_t *) ML_allocate( tsz2*sizeof(int_t) );
      lugrid_tiles = (ML_Lugrid *) ML_allocate( ntile*sizeof(ML_Lugrid) );
      k = 0;
      for( g=0 ; g < ntile; g++)
      {
         if( g < stile )
         {
            tsz2 = (tsz-1)*(tsz-1);
            nprow = tsz-1;
            npcol = tsz-1;
         }
         else if( g < stile+mtile)
         {
            tsz2 = tsz*(tsz-1);
            nprow = tsz;
            npcol = tsz-1;
         }
         else
         {
            tsz2 = tsz*tsz;
            nprow = tsz;
            npcol = tsz;
         }
         for( l=0; l<tsz2; l++)
         {
            usermap[l] = l+k;
            if( iam == l+k ) mygroup = g;
         }
         k = k + tsz2;
         /* in-lining
          * superlu_gridmap( MPI_COMM_WORLD, 
          * nprow, npcol, usermap, nprow, &((lugrid_tiles[g]).grid));
          */
         (lugrid_tiles[g]).grid.nprow = nprow;
         (lugrid_tiles[g]).grid.npcol = npcol;
         MPI_Comm_group( MPI_COMM_WORLD, &mpi_base_group );
         MPI_Group_incl( mpi_base_group, tsz2, usermap, &superlu_grp );
         MPI_Comm_create(MPI_COMM_WORLD,superlu_grp,&(lugrid_tiles[g].grid.comm));
         if ( lugrid_tiles[g].grid.comm == MPI_COMM_NULL ) 
         {
            lugrid_tiles[g].grid.comm = MPI_COMM_WORLD;
            lugrid_tiles[g].grid.iam  = iam;
         } 
         else 
         {   /* mygroup=g and iam%mtile=npcol*myrow + mycol */
            MPI_Comm_rank(lugrid_tiles[g].grid.comm,&(lugrid_tiles[g].grid.iam));
            myrow = lugrid_tiles[g].grid.iam / npcol;
            mycol = lugrid_tiles[g].grid.iam % npcol;
            MPI_Comm_split(lugrid_tiles[g].grid.comm, 
                           myrow, mycol, &(lugrid_tiles[g].grid.rscp.comm));
            MPI_Comm_split(lugrid_tiles[g].grid.comm, 
                           mycol, myrow, &(lugrid_tiles[g].grid.cscp.comm));
            lugrid_tiles[g].grid.rscp.Np  = npcol;
            lugrid_tiles[g].grid.rscp.Iam = mycol;
            lugrid_tiles[g].grid.cscp.Np  = nprow;
            lugrid_tiles[g].grid.cscp.Iam = myrow;
         }
      } /* end for group g */
      if( nprocs != k )
      {
         printf("Error nprocs %d  k %d \n", nprocs, k);
         exit(-1);
      }
      ML_free(usermap);
      solver->ML_subgroup = mygroup;
      /*
       * Fact = DOFACT Trans = NOTRANS Equil = EQUI RowPerm = LargeDiag
       * ColPerm = COLAMD ReplaceTinyPivot = REPLACE IterRefine = DOUBLE
       */
      set_default_options(optionsptr);
      optionsptr->Equil = NEQU;
      optionsptr->IterRefine = NOREFINE;
      /*
       * Future possiblities to experiment with include
       * optionsptr->RowPerm = NOROWPERM;
       * optionsptr->ColPerm = COLAMD;   (default)
       * optionsptr->ColPerm = MMD_AT_PLUS_A;
       * optionsptr->ColPerm = NATURAL;
       * ... and Equil
       */
      ScalePermstruct = ( ScalePermstruct_t *) ML_allocate( sizeof( ScalePermstruct_t));
      ScalePermstructInit(n, n, ScalePermstruct);
      LUstruct = ( LUstruct_t *) ML_allocate( sizeof( LUstruct_t) );
      LUstructInit(n, n, LUstruct);
      solver->PERMspl = (void *) ScalePermstruct;
      solver->LUspl = (void *) LUstruct;
      solver->gridtiles = lugrid_tiles;
      /* rst: mygrid is a pointer to a structure, not a structure.
       *  mygrid = ( gridinfo_t *) ML_allocate( sizeof( gridinfo_t) );
       */
   } 
   else 
   {
      /* Indicate that the factored form of A is supplied. */
      /* Reset options */
      optionsptr->Fact = FACTORED;
      optionsptr->Trans = NOTRANS;
      optionsptr->Equil = NEQU;
      optionsptr->RowPerm = MY_PERMR;
      optionsptr->ColPerm = MY_PERMC;
      optionsptr->ReplaceTinyPivot = REPLACE;
      optionsptr->IterRefine = NOREFINE;
      lugrid_tiles    = solver->gridtiles;
      mygroup = (int ) solver->ML_subgroup;
   }
   mygrid = &((lugrid_tiles[mygroup]).grid);

   /* ------------------------------------------------------------- */
   /* gather from all processors the complete right hand side       */
   /* ------------------------------------------------------------- */

   nrhs = 1;
   ML_memory_alloc((void**) &local_rhs, n*sizeof(double),"LU1" );
   for ( i = 0; i < N_local; i++ ) local_rhs[i] = rhs[i];
   i = N_local;
   ML_Comm_GappendDouble((ML_Comm *) mlcomm, local_rhs, &i, n);

   /* ------------------------------------------------------------- */
   /* perform LU decomposition and then solve                       */
   /* ------------------------------------------------------------- */

   info = flag;
   PStatInit(&stat);
   pdgssvx_ABglobal(optionsptr, A, ScalePermstruct, local_rhs, n, 
                    nrhs, mygrid, LUstruct, &berr, &stat, &info);

/*
if( iam == 0 ){
heap_info(&fragments, &total_free, &largest_free, &total_used);
printf("memory usage: fragments %d free: total %d, largest %d, total_used %d\n",
          fragments, total_free, largest_free, total_used);
}
*/

   if ( flag == 0 ) 
   {
     if ( A != NULL ) 
     {
/*
        ML_memory_free((void*)&(((NCformat *) A->Store)->rowind));
        ML_memory_free((void*)&(((NCformat *) A->Store)->colptr));
        ML_memory_free((void*)&(((NCformat *) A->Store)->nzval));
        SUPERLU_FREE( A->Store );
*/
        /* to satisfy pdgssvx_ABglobal argument check, postpone
         * ML_memory_free((void**) &A);
         */
      }
   }

   solver->reuse_flag = 1;
   PStatFree(&stat);
   if( info != 0 )
   {
      if( iam == 0 )printf("Error: ml_superlu    info = %d\n",info);
      return(-1);
   }
   /* ------------------------------------------------------------- */
   /* extract the local solution sub-vector and then clean up       */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < N_local; i++ ) x[i] = local_rhs[i+offset];
   ML_memory_free( (void **) &local_rhs );
#else
   printf("ML_SuperLU_Solve : SuperLU not used.\n");
   ML_avoid_unused_param( (void *) vsolver);
   ML_avoid_unused_param( (void *) &ilen);
   ML_avoid_unused_param( (void *) x);
   ML_avoid_unused_param( (void *) &olen);
   ML_avoid_unused_param( (void *) rhs);
#endif
#endif
   return 0;
}

/* ******************************************************************** */
/* This subroutine calls the SuperLU subroutine to solve a given        */
/* subproblem where the matrix and right hand side are both residing    */
/* in the local processor (domain decomposition)                        */
/* ******************************************************************** */

int ML_SuperLU_SolveLocal(void *vsolver, double *x, double *rhs)
{
#ifdef ML_SUPERLU2

   SuperMatrix *A, *L, *U, B;
   int         i, n, info, flag, *perm_r, *perm_c, permc_spec;
   ML_Solver   *solver;

   /* ------------------------------------------------------------- */
   /* fetch the sparse matrix and other parameters                  */
   /* ------------------------------------------------------------- */

   solver  = (ML_Solver *) vsolver;
   A       = (SuperMatrix *) solver->Mat1;
   n       = (int) solver->dble_params1[0];
   flag    = solver->reuse_flag;
   if ( flag != 0 ) 
   {
      for ( i = 0; i < n; i++ ) x[i] = rhs[i];
   }

   /* ------------------------------------------------------------- */
   /* if factorization has not been done, allocate space for it     */
   /* ------------------------------------------------------------- */

   if ( flag == 0 ) 
   {
      ML_memory_alloc((void**) &perm_c, n * sizeof(int), "LU6" );
      ML_memory_alloc((void**) &perm_r, n * sizeof(int), "LU7" );
      solver->int_params1 = perm_r;
      solver->int_params2 = perm_c;
      permc_spec = 0;
      get_perm_c(permc_spec, A, perm_c);
      ML_memory_alloc((void**) &L, sizeof(SuperMatrix), "LU8" );
      ML_memory_alloc((void**) &U, sizeof(SuperMatrix), "LU9" );
      solver->Mat2 = (void *) L;
      solver->Mat3 = (void *) U;
   } 
   else 
   {
      perm_r = (int *) solver->int_params1;
      perm_c = (int *) solver->int_params2;
      L = (SuperMatrix *) solver->Mat2;
      U = (SuperMatrix *) solver->Mat3;
   }

   /* ------------------------------------------------------------- */
   /* create right hand side vector conforming to SuperLU format    */
   /* ------------------------------------------------------------- */

   dCreate_Dense_Matrix(&B, n, 1, x, n, DN, _D, GE);
   info = flag;
   dgssv(A, perm_c, perm_r, L, U, &B, &info);

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   Destroy_SuperMatrix_Store(&B);
   solver->reuse_flag = 1;
#else
   printf("ML_SuperLU_SolveLocal : SuperLU not used.\n");
   ML_avoid_unused_param( (void *) vsolver);
   ML_avoid_unused_param( (void *) x);
   ML_avoid_unused_param( (void *) rhs);
#endif
   return 0;
}

/* ************************************************************************* */
/* clean up                                                                  */
/* ------------------------------------------------------------------------- */

int ML_CSolve_Clean_SuperLU( void *vsolver, ML_CSolveFunc *func)
{
   ML_Solver   *solver;

#ifdef SUPERLU
   SuperMatrix *Amat;

   solver = (ML_Solver *) vsolver;
   solver->reuse_flag = -999;
   func->internal( vsolver, 0, NULL, 0, NULL);

   Amat = (SuperMatrix*) solver->Mat1;
   if (Amat != NULL )
   {
      SUPERLU_FREE( ((NRformat *) Amat->Store)->colind);
      SUPERLU_FREE( ((NRformat *) Amat->Store)->rowptr);
      SUPERLU_FREE( ((NRformat *) Amat->Store)->nzval);
      SUPERLU_FREE( Amat->Store );
      ML_memory_free(  (void**) &(solver->Mat1) );
      solver->Mat1 = NULL;
   }
#else
#ifdef DSUPERLU
   SuperMatrix *Amat;

   solver = (ML_Solver *) vsolver;
   solver->reuse_flag = -999;
   func->internal( vsolver, 0, NULL, 0, NULL);
   Amat = (SuperMatrix*) solver->Mat1;
   if (Amat != NULL)
   {
      Destroy_CompCol_Matrix(Amat);
      ML_memory_free((void**) &Amat);
   }
   solver->Mat1 = NULL;

#else
   solver = (ML_Solver *) vsolver;
   solver->reuse_flag = -999;
   func->internal( vsolver, 0, NULL, 0, NULL);
#endif
#endif
   ML_Solver_Destroy( &solver );
   return 0;
}

/* ******************************************************************** */
/* An array of processors decomposes into tiles of three sizes:         */
/* (n-1)^2 , n (n-1) or n^2.  * For any natural numbers, p and n,       */
/* there exist integers i,j,k such that 0<=i,j<n and                    */
/* p = i*(n-1)^2 + j*n*(n-1) + k*n^2 = (i+j+k)*(n^2) -n*(2*i+j)+i       */
/* If p >= (2*n-1)*(n-1)^2 , then k >= 0                                */
/* ******************************************************************** */

void ML_SuperLU_Set_Tile( int p, int* n, int* i, int* j, int* k)
{
   int l,q,r,s;
   double quotient, cuberoot; 
   if( p < 12) 
   {
      *n = 2;
      *i = p;
      *j = 0;
      *k = 0;
   }
   else if( p < 54) 
   {
      /* p = l + 4s + 12 q */
      l = p % 4;
      r = (p-l)/4;
      s = r % 3;
      q = (r-s)/3;
      if( l == 0 )
      {
         *n = 2;
         *i = 0;
         *j = 0;
         *k = s + 3*q;
      }
      else if( l == 1 )
      {
         *n = 3;
         *i = s + 3*q -2;
         *j = 0;
         *k = 1; 
      }
      else if( l == 2 )
      {
         *n = 3;
         *i = s + 3*q -1;
         *j = 1;
         *k = 0; 
      }
      else
      { 
         *n = 3;
         *i = s + 3*q -3;
         *j = 1;
         *k = 1; 
      }
   }
   else
   {
      quotient = (double) p / 2;
      cuberoot = (double) 1 / 3;
      *n = (int) floor(pow(quotient,cuberoot));
      *i = p % (*n);
      q = (p - *i )/ *n;
      r = (q + 2 * *i ) % (*n);
      if( r > 0 )
         *j = *n - r;
      else
         *j = 0;
      *k = ( (q + 2 * *i + *j ) / (*n) ) - *i - *j;
   }
} /* end of ML_SuperLU_Set_Tile */


/*****************************************************************************/
/* clean up                                                                  */
/* ------------------------------------------------------------------------- */

int ML_Clean_CSolveSuperLU( void *vsolver, ML_CSolveFunc *func) 
{
   ML_Solver   *solver;

#ifdef SUPERLU
   SuperMatrix *Amat;

   solver = (ML_Solver *) vsolver;
   solver->reuse_flag = -999;
   func->internal( vsolver, 0, NULL, 0, NULL);

   Amat = (SuperMatrix*) solver->Mat1;
   if (Amat != NULL ) {
      SUPERLU_FREE( ((NRformat *) Amat->Store)->colind);
      SUPERLU_FREE( ((NRformat *) Amat->Store)->rowptr);
      SUPERLU_FREE( ((NRformat *) Amat->Store)->nzval);
      SUPERLU_FREE( Amat->Store );
      ML_memory_free(  (void**) &(solver->Mat1) );
      solver->Mat1 = NULL;
   }
#else
#ifdef DSUPERLU
   SuperMatrix *Amat;

   solver = (ML_Solver *) vsolver;
   solver->reuse_flag = -999;
   func->internal( vsolver, 0, NULL, 0, NULL);
   Amat = (SuperMatrix*) solver->Mat1;
   if (Amat != NULL) {
      Destroy_CompCol_Matrix(Amat);
      ML_memory_free((void**) &Amat);
   }
   solver->Mat1 = NULL;

#else
   solver = (ML_Solver *) vsolver;
   solver->reuse_flag = -999;
   func->internal( vsolver, 0, NULL, 0, NULL);
#endif
#endif
   ML_Solver_Destroy( &solver );
   return 0;
}

/*****************************************************************************/
/* Generate a coarse grid matrix suitable for solution with SuperLU          */
/* ------------------------------------------------------------------------- */

int ML_Gen_CoarseSolverSuperLU(ML *ml_handle, int level)
{
#ifdef SUPERLU
   int            i, j, *mat_ia, *mat_ja, nrows, nnz, offset, N_local;
   int            reuse, coarsest_level, flag, space, *cols, nz_ptr;
   int            getrow_flag, osize, *row_ptr, length, zero_flag;
   double         *mat_val, *vals, dsize, di;
   void           *data;
   ML_1Level      *sl;
   ML_Solver      *solver;
   ML_Operator    *op;
   SuperMatrix    *A;
   ML_Matrix_DCSR *csr_mat, *csr2_mat;
struct ML_CSR_MSRdata *temp_ptr;
ML *subml;
int nblocks = 1, *block_list, old_upper = 0, count, newptr, me, nnzs;
#ifdef ML_TIMING
   double t0;

   t0 = GetClock();
#endif

   /* ----------------------------------------------------------------- */
   /* extract local matrix using getrow function and store it into a    */
   /* CSR data object                                                   */
   /* ----------------------------------------------------------------- */

   if ( level < 0 || level >= ml_handle->ML_num_levels ) {
      printf("ML_Gen_CoarseSolverSuperLU error : invalid level number.\n");
      exit(-1);
   }
   op      = (ML_Operator *) &ml_handle->Amat[level];
   data    = op->data;
   osize   = op->outvec_leng;
   if (op->invec_leng < 0) {
      nblocks = -op->invec_leng;
      op->invec_leng = osize;
   }
   row_ptr = (int *) ML_allocate(sizeof(int)*(osize+1));
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
      cols    = (int    *) ML_allocate(sizeof(int)*space);
      vals    = (double *) ML_allocate(sizeof(double)*space);

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
   csr_mat = (ML_Matrix_DCSR *) ML_allocate(sizeof(ML_Matrix_DCSR));
   csr_mat->mat_n  = osize;
   csr_mat->mat_ja = cols;
   csr_mat->mat_a  = vals;
   csr_mat->mat_ia = row_ptr;
   csr_mat->comminfo = op->getrow->pre_comm;

   /* ----------------------------------------------------------------- */
   /* form a global matrix                                              */
   /* ----------------------------------------------------------------- */

   csr2_mat = (ML_Matrix_DCSR *) ML_allocate(sizeof(ML_Matrix_DCSR));
   ML_Gen_Amatrix_Global( csr_mat, csr2_mat, ml_handle->comm, &offset);
   ML_free(row_ptr);
   ML_free(cols);
   ML_free(vals);
   ML_free(csr_mat);

   /* Throw away some information to make it cheaper for LU. We do this   */ 
   /* by using metis to generate some blocks and factor the block matrix. */
   if (nblocks > 1) { 
      mat_ia  = csr2_mat->mat_ia;
      mat_ja  = csr2_mat->mat_ja;
      mat_val = csr2_mat->mat_a;
      nrows   = csr2_mat->mat_n;
      temp_ptr =(struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
      temp_ptr->rowptr = mat_ia;
      temp_ptr->columns= mat_ja;
      temp_ptr->values = mat_val;
      ML_Create(&subml, 1);
      ML_Init_Amatrix(subml, 0, nrows, nrows, (void *) temp_ptr);
      MLnew_Set_Amatrix_Matvec(subml, 0, CSR_matvec);
      ML_CommInfoOP_Set_neighbors(&(subml->Amat[0].getrow->pre_comm), 0,
                               NULL, ML_OVERWRITE, NULL, 0);
      ML_Operator_Set_Getrow(&(subml->Amat[0]), ML_INTERNAL, 
                             subml->Amat[0].outvec_leng, CSR_getrow);
      ML_Gen_Blocks_Metis(subml, 0, &nblocks, &block_list);
      ML_Destroy(&subml);
      ML_free(temp_ptr);
      for (i = 0; i < nrows; i++) {
         me = block_list[i];
         for (j = mat_ia[i]; j < mat_ia[i+1]; j++) {
            if ( block_list[mat_ja[j]] != me) {mat_ja[j] = -1; }
         }
      }
      ML_free(block_list);

      if (nrows > 0) old_upper = mat_ia[0];
      nnzs = mat_ia[nrows];
      for (i = 0; i < nrows; i++) {
	count = 0;
        for (j = old_upper; j < mat_ia[i+1]; j++) {
           if ( mat_ja[j] != -1) count++;
        }
        old_upper = mat_ia[i+1];
        mat_ia[i+1] = mat_ia[i] + count;
      }

      newptr = 0;
      for (i = 0; i < nnzs; i++) {
         if ( mat_ja[i] != -1) {
            mat_ja[newptr] = mat_ja[i];
            mat_val[newptr++] = mat_val[i];
         }
      }
   }


   /* ----------------------------------------------------------------- */
   /* set SuperLU as solver                                             */
   /* ----------------------------------------------------------------- */

   coarsest_level = level;
   sl = &(ml_handle->SingleLevel[coarsest_level]);
   if ( sl->csolve->func->internal == ML_SuperLU_Solve ) reuse = 1;
   else
   {
      reuse = 0;
      sl->csolve->func->internal = ML_SuperLU_Solve;
      sl->csolve->func->ML_id = ML_INTERNAL;
      ML_CSolve_Set_Label( sl->csolve, "SuperLU");

   }

   /* ----------------------------------------------------------------- */
   /* free up previous storage                                          */
   /* ----------------------------------------------------------------- */

   if ( sl->csolve->data != NULL )
   {
      solver = (ML_Solver *) sl->csolve->data;
      if ( reuse == 1 )
      {
         /* Charles look at these  */
        /* if (solver->int_params1 != NULL)
         {
            ML_memory_free( (void**) &(solver->int_params1) );
            solver->int_params1 = NULL;
         }
         if (solver->int_params2 != NULL)
         {
            ML_memory_free( (void**) &(solver->int_params2) );
            solver->int_params2 = NULL;
         }*/
         if ( solver->dble_params1 != NULL )
         {
            ML_memory_free(  (void**) &(solver->dble_params1) );
            solver->dble_params1 = NULL;
         }
         solver->reuse_flag = -999;
         ML_SuperLU_Solve((void*)solver, 0, NULL, 0, NULL);
         solver->reuse_flag = 0;
         /* Charles look at these  */
         /* if (solver->Mat1 != NULL )
         {
            Destroy_CompRow_Matrix(solver->Mat1);
            ML_memory_free(  (void**) &(solver->Mat1) ); 
            solver->Mat1 = NULL;
         }
         if (solver->Mat2 != NULL )
         {
            Destroy_SuperNode_Matrix(solver->Mat2);
            ML_memory_free(  (void**) &(solver->Mat2) ); 
            solver->Mat2 = NULL;
         }
         if (solver->Mat3 != NULL )
         {
            Destroy_CompCol_Matrix(solver->Mat3);
            ML_memory_free(  (void**) &(solver->Mat3) ); 
            solver->Mat3 = NULL;
         }*/
      }
      ML_memory_free(  (void**) &(solver) );
   }

   /* ----------------------------------------------------------------- */
   /* create new context                                                */
   /* ----------------------------------------------------------------- */

   ML_Solver_Create( &solver );
   sl->csolve->data = (void *) solver;
   solver->reuse_flag = 0;
   solver->void_params1 = (void *) ml_handle->comm;
   ML_memory_alloc( (void **) &vals, 3 * sizeof(double), "KLI" );
   N_local = osize;
   vals[0]  = (double) N_local;
   vals[1]  = (double) offset;
   vals[2]  = (double) csr2_mat->mat_n;
   solver->dble_params1 = (double *) vals;

   /* ----------------------------------------------------------------- */
   /* form SuperLU type matrix                                          */
   /* ----------------------------------------------------------------- */

   mat_ia  = csr2_mat->mat_ia;
   mat_ja  = csr2_mat->mat_ja;
   mat_val = csr2_mat->mat_a;
   nrows   = csr2_mat->mat_n;
   nnz     = mat_ia[nrows];
   ML_memory_alloc( (void **) &A, sizeof(SuperMatrix), "KLJ" );
   dCreate_CompRow_Matrix(A,nrows,nrows,nnz,mat_val,mat_ja,mat_ia,NR,_D,GE);
   solver->Mat1 = (void *) A;
   /* Charles look at these */
   /* solver->Mat1 = NULL;
   SUPERLU_FREE(A->Store);
   ML_memory_free( (void **) &A );
   ML_memory_free( (void **) &mat_ia );
   ML_memory_free( (void **) &mat_ja );
   ML_memory_free( (void **) &mat_val ); */
   ML_free(csr2_mat);
#ifdef ML_TIMING
   sl->csolve->build_time = GetClock() - t0;
   ml_handle->timing->total_build_time += sl->csolve->build_time;
#endif

#else
#ifdef DSUPERLU
   int               i, offset, N_local;
   int               reuse, coarsest_level, flag, space, *cols, nz_ptr;
   int               getrow_flag, osize, *row_ptr, length;
   int               j, k, k1, k2, next,*ia, *ja;
   int_t             *mat_ia, *mat_ja, nrows, nnz;
   double            *mat_val, *vals, dsize, di, *aa;
   void              *data;
   ML_1Level         *sl;
   ML_Solver         *solver;
   ML_Operator       *op;
   SuperMatrix       *A;
   ML_Matrix_DCSR    *csr_mat, *csr2_mat;
   struct ML_CSR_MSRdata *temp_ptr;
int nblocks = 1, *block_list, old_upper = 0, count, newptr, me, nnzs;
   ML *subml;

   /* ----------------------------------------------------------------- */
   /* extract local matrix using getrow function and store it into a    */
   /* CSR data object                                                   */
   /* ----------------------------------------------------------------- */

   if ( level < 0 || level >= ml_handle->ML_num_levels ) 
   {
      printf("ML_Gen_CoarseSolverSuperLU error : invalid level number.\n");
      exit(-1);
   }
   op      = (ML_Operator *) &ml_handle->Amat[level];
   data    = op->data;
   osize   = op->outvec_leng;
   if (op->invec_leng < 0) 
   {
      nblocks = -op->invec_leng;
      op->invec_leng = osize;
   }
   row_ptr = (int *) ML_allocate(sizeof(int)*(osize+1));
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
      cols    = (int    *) ML_allocate(sizeof(int)*space);
      vals    = (double *) ML_allocate(sizeof(double)*space);

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
   csr_mat = (ML_Matrix_DCSR *) ML_allocate(sizeof(ML_Matrix_DCSR));
   csr_mat->mat_n  = osize;
   csr_mat->mat_ja = cols;
   csr_mat->mat_a  = vals;
   csr_mat->mat_ia = row_ptr;
   csr_mat->comminfo = op->getrow->pre_comm;

   /* ----------------------------------------------------------------- */
   /* form a global matrix                                              */
   /* SuperLU_Dist support has been curtailed, particularly for NR      */
   /* mat := csr2_mat (in column format) = csr2_mat transpose           */
   /* ----------------------------------------------------------------- */

   csr2_mat = (ML_Matrix_DCSR *) ML_allocate(sizeof(ML_Matrix_DCSR));
   ML_Gen_Amatrix_Global( csr_mat, csr2_mat, ml_handle->comm, &offset);
   ML_free(cols);
   ML_free(vals);
   ML_free(row_ptr);
   ML_free(csr_mat);

   /* Throw away some information to make it cheaper for LU. We do this   */ 
   /* by using metis to generate some blocks and factor the block matrix. */
   if (nblocks > 1) { 
      mat_ia  = csr2_mat->mat_ia;
      mat_ja  = csr2_mat->mat_ja;
      mat_val = csr2_mat->mat_a;
      nrows   = csr2_mat->mat_n;
      temp_ptr =(struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
      temp_ptr->rowptr = mat_ia;
      temp_ptr->columns= mat_ja;
      temp_ptr->values = mat_val;
      ML_Create(&subml, 1);
      ML_Init_Amatrix(subml, 0, nrows, nrows, (void *) temp_ptr);
      ML_CommInfoOP_Set_neighbors(&(subml->Amat[0].getrow->pre_comm), 0,
                               NULL, ML_OVERWRITE, NULL, 0);
      ML_Operator_Set_Getrow(&(subml->Amat[0]), ML_INTERNAL, 
                             subml->Amat[0].outvec_leng, CSR_getrow);

      MLnew_Set_Amatrix_Matvec(subml, 0, CSR_matvec);
      ML_Gen_Blocks_Metis(subml, 0, &nblocks, &block_list);
      ML_Destroy(&subml);
      ML_free(temp_ptr);
      for (i = 0; i < nrows; i++) {
         me = block_list[i];
         for (j = mat_ia[i]; j < mat_ia[i+1]; j++) {
            if ( block_list[mat_ja[j]] != me) {mat_ja[j] = -1; }
         }
      }
      ML_free(block_list);

      if (nrows > 0) old_upper = mat_ia[0];
      nnzs = mat_ia[nrows];
      for (i = 0; i < nrows; i++) {
	count = 0;
        for (j = old_upper; j < mat_ia[i+1]; j++) {
           if ( mat_ja[j] != -1) count++;
        }
        old_upper = mat_ia[i+1];
        mat_ia[i+1] = mat_ia[i] + count;
      }

      newptr = 0;
      for (i = 0; i < nnzs; i++) {
         if ( mat_ja[i] != -1) {
            mat_ja[newptr] = mat_ja[i];
            mat_val[newptr++] = mat_val[i];
         }
      }
   }

   /*
    * if (global_comm->ML_mypid == 0) {
    *  for (i = 0; i <= csr2_mat->mat_n; i++)
    *   printf("row_ptr(%d) = %d\n",i,csr2_mat->mat_ia[i]);
    *  for (i = 0; i < csr2_mat->mat_ia[csr2_mat->mat_n]; i++)
    *   printf("(%d,   %d,%e)\n",i,csr2_mat->mat_ja[i],csr2_mat->mat_a[i]);
    * }
    */
   nrows   = csr2_mat->mat_n;
   nnz     = csr2_mat->mat_ia[nrows];
   ia      = csr2_mat->mat_ia;
   ja      = csr2_mat->mat_ja;
   aa      = csr2_mat->mat_a;
   ML_memory_alloc( (void **) &mat_val, nnz*sizeof(double), "cat" );
   ML_memory_alloc( (void **) &mat_ja, nnz*sizeof(int), "jct" );
   ML_memory_alloc( (void **) &mat_ia, (nrows+1)*sizeof(int), "ict" );
   for(i=0;i<=nrows;i++) mat_ia[i] = 0;
   for(i=0;i<nrows;i++){
     k1 = ia[i];
     k2 = ia[i+1];
     for(k=k1;k<k2;k++){
       j = ja[k]+1;
       ++mat_ia[j];
     }
   }
   for(i=0;i<nrows;i++)mat_ia[i+1] = mat_ia[i] + mat_ia[i+1];
   for(i=0;i<nrows;i++){
     k1 = ia[i];
     k2 = ia[i+1];
     for(k=k1;k<k2;k++){
       j = ja[k];
       next = mat_ia[j];
       mat_ia[j] = next+1;
       mat_ja[next] = i;
       mat_val[next] = aa[k];
     }
   }
   for(i=nrows-1;i>=0;i--)mat_ia[i+1] = mat_ia[i];
   mat_ia[0] = 0;
   ML_memory_free(  (void**) &(csr2_mat->mat_ia) );
   ML_memory_free(  (void**) &(csr2_mat->mat_ja) );
   ML_memory_free(  (void**) &(csr2_mat->mat_a) );
   csr2_mat->mat_ia = mat_ia;
   csr2_mat->mat_ja = mat_ja;
   csr2_mat->mat_a  = mat_val;

   /* ----------------------------------------------------------------- */
   /* set SuperLU as solver                                             */
   /* ----------------------------------------------------------------- */

   coarsest_level = level;
   sl = &(ml_handle->SingleLevel[coarsest_level]);
   if ( sl->csolve->func->internal == ML_SuperLU_Solve ) reuse = 1;
   else
   {
      reuse = 0;
      sl->csolve->func->internal = ML_SuperLU_Solve;
      sl->csolve->func->ML_id = ML_INTERNAL;
      ML_CSolve_Set_Label( sl->csolve, "Dist. SuperLU");
   }

   /* ----------------------------------------------------------------- */
   /* free up previous storage                                          */
   /* ----------------------------------------------------------------- */

   if ( sl->csolve->data != NULL )
   {
      solver = (ML_Solver *) sl->csolve->data;
      if ( reuse == 1 )
      {
         if (solver->int_params1 != NULL)
         {
            ML_memory_free( (void**) &(solver->int_params1) );
            solver->int_params1 = NULL;
         }
         if (solver->int_params2 != NULL)
         {
            ML_memory_free( (void**) &(solver->int_params2) );
            solver->int_params2 = NULL;
         }
         if ( solver->dble_params1 != NULL )
         {
            ML_memory_free(  (void**) &(solver->dble_params1) );
            solver->dble_params1 = NULL;
         }
         if (solver->Mat1 != NULL )
         {
            Destroy_CompCol_Matrix(solver->Mat1);
            ML_memory_free(  (void**) &(solver->Mat1) );
            solver->Mat1 = NULL;
         }
         if (solver->Mat2 != NULL )
         {
            Destroy_SuperNode_Matrix(solver->Mat2);
            ML_memory_free(  (void**) &(solver->Mat2) );
            solver->Mat2 = NULL;
         }
         if (solver->Mat3 != NULL )
         {
            Destroy_CompCol_Matrix(solver->Mat3);
            ML_memory_free(  (void**) &(solver->Mat3) );
            solver->Mat3 = NULL;
         }
      }
      ML_memory_free(  (void**) &(solver) );
   }
   /* ----------------------------------------------------------------- */
   /* create new context                                                */
   /* ----------------------------------------------------------------- */

   ML_Solver_Create( &solver );
   sl->csolve->data = (void *) solver;
   solver->reuse_flag = 0;
   solver->void_params1 = (void *) ml_handle->comm;
   ML_memory_alloc( (void **) &vals, 3 * sizeof(double), "KLI" );
   N_local = osize;
   vals[0]  = (double) N_local;
   vals[1]  = (double) offset;
   vals[2]  = (double) csr2_mat->mat_n;
   solver->dble_params1 = (double *) vals;

   /* ----------------------------------------------------------------- */
   /* form SuperLU type matrix                                          */
   /* ----------------------------------------------------------------- */

   ML_memory_alloc( (void **) &A, sizeof(SuperMatrix), "KLJ" );
   dCreate_CompCol_Matrix(A,nrows,nrows,nnz,mat_val,mat_ja,mat_ia,NC,_D,GE);
   solver->Mat1 = (void *) A;
   ML_free(csr2_mat);
#else
   printf("ML : SuperLU not linked.\n");
   ML_avoid_unused_param( (void *) ml_handle);
   ML_avoid_unused_param( (void *) &level);
#endif
#endif
   return 0;
}

/* ************************************************************************* */
/* variable block additive Schwarz                                           */
/* ------------------------------------------------------------------------- */

int ML_Smoother_VBlockAdditiveSchwarz(ML_Smoother *sm, int inlen, double x[],
                                      int outlen, double rhs[])
{
#ifdef SUPERLU
   int                i, j, m, k, extNrows, nblocks, length, *indptr, ntimes;
   int                *blk_size, **blk_indices, max_blk_size;
   int                **aux_bmat_ia, **aux_bmat_ja;
   double             *dbuffer, **aux_bmat_aa, *rhsbuf, *solbuf, *xbuffer = NULL;
   ML_Comm            *comm;
   ML_Operator        *Amat;
   ML_Smoother        *smooth_ptr;
   ML_Sm_Schwarz_Data *dataptr;
   ML_CommInfoOP      *getrow_comm;
   int                info, *perm_r, *perm_c, *etree, panel_size, lwork;
   double             *R, *C, *ferr, *berr, rpg, rcond, dtemp;
   char               fact[1], equed[1], trans[1], refact[1];
   void               *work=NULL;
   SuperMatrix        *A, *L, *U, B, X;
   factor_param_t     iparam;
   mem_usage_t        mem_usage;

   /* --------------------------------------------------------- */
   /* fetch parameters and check                                */
   /* --------------------------------------------------------- */

   smooth_ptr = (ML_Smoother *) sm;
   comm       = smooth_ptr->my_level->comm;
   Amat       = smooth_ptr->my_level->Amat;
   dataptr    = (ML_Sm_Schwarz_Data *) smooth_ptr->smoother->data;
   ntimes     = smooth_ptr->ntimes;

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_Smoother_AdditiveSchwarz): Need getrow()\n");
   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for AdditiveSchwarz\n");
   if ( dataptr == NULL ) 
      pr_error("Error(AdditiveSchwarz): Need dataptr\n");

   getrow_comm = Amat->getrow->pre_comm;
   extNrows    = dataptr->Nrows;
   nblocks     = dataptr->nblocks;
   blk_indices = dataptr->blk_indices;
   blk_size    = dataptr->blk_size;
   aux_bmat_ia = dataptr->aux_bmat_ia;
   aux_bmat_ja = dataptr->aux_bmat_ja;
   aux_bmat_aa = dataptr->aux_bmat_aa;
   max_blk_size = 0;
   for ( i = 0; i < nblocks; i++ ) 
      max_blk_size = (blk_size[i] > max_blk_size) ? blk_size[i] : max_blk_size;

   /* --------------------------------------------------------- */
   /* communicate the rhs and put into dbuffer                  */
   /* --------------------------------------------------------- */

   dbuffer = (double *) ML_allocate(extNrows * sizeof(double));
   for ( i = 0; i < outlen; i++ ) dbuffer[i] = rhs[i];
   for ( i = 0; i < inlen;  i++ ) x[i] = 0.0;

   if (extNrows > outlen && getrow_comm != NULL)
      ML_exchange_bdry(dbuffer,getrow_comm,inlen,comm,ML_OVERWRITE,NULL);

   /* --------------------------------------------------------- */
   /* set up for SuperLU solves                                 */
   /* --------------------------------------------------------- */

   rhsbuf = (double *) ML_allocate(max_blk_size * sizeof(double));
   solbuf = (double *) ML_allocate(max_blk_size * sizeof(double));
   panel_size               = sp_ienv(1);
   iparam.panel_size        = panel_size;
   iparam.relax             = sp_ienv(2);
   iparam.diag_pivot_thresh = 1.0;
   iparam.drop_tol          = -1;
   lwork                    = 0;
   *fact                    = 'F';
   *equed                   = 'N';
   *trans                   = 'N';
   *refact                  = 'N';
   R     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   C     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   ferr  = (double *) SUPERLU_MALLOC(sizeof(double));
   berr  = (double *) SUPERLU_MALLOC(sizeof(double));
   etree = (int *) ML_allocate( max_blk_size * sizeof(int) );

   /* --------------------------------------------------------- */
   /* the first pass                                            */
   /* --------------------------------------------------------- */

   for ( i = 0; i < nblocks; i++ )
   {
      indptr = blk_indices[i];
      length = blk_size[i];
      for ( j = 0; j < length; j++ ) rhsbuf[j] = dbuffer[indptr[j]];
      A = dataptr->slu_Amat[i];
      L = dataptr->slu_Lmat[i];
      U = dataptr->slu_Umat[i];
      perm_c = dataptr->perm_c[i];
      perm_r = dataptr->perm_r[i];
      dCreate_Dense_Matrix(&B, length, 1, rhsbuf, length, DN, _D, GE);
      dCreate_Dense_Matrix(&X, length, 1, solbuf, length, DN, _D, GE);
      dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
             equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
             ferr, berr, &mem_usage, &info);
      for ( j = 0; j < length; j++ ) 
         /*if ( indptr[j] < inlen ) x[indptr[j]] += solbuf[j];*/
if ( indptr[j] < inlen ) x[indptr[j]] = solbuf[j];
      Destroy_SuperMatrix_Store(&B);
      Destroy_SuperMatrix_Store(&X);
   }

   if (ntimes > 1) xbuffer = (double *) ML_allocate(extNrows * sizeof(double));

   for ( m = 1; m < ntimes; m++ )
   {
      for ( i = 0; i < inlen; i++ ) xbuffer[i] = x[i];
      if (extNrows > outlen && getrow_comm != NULL)
         ML_exchange_bdry(xbuffer,getrow_comm,inlen,comm,ML_OVERWRITE,NULL);

      for ( i = 0; i < nblocks; i++ )
      {
         indptr = blk_indices[i];
         length = blk_size[i];
         for ( j = 0; j < length; j++ )
         {
            dtemp = dbuffer[indptr[j]];
            for ( k = aux_bmat_ia[i][j]; k < aux_bmat_ia[i][j+1]; k++ )
               dtemp -= (aux_bmat_aa[i][k] * xbuffer[aux_bmat_ja[i][k]]); 
            rhsbuf[j] = dtemp;
         }
         A = dataptr->slu_Amat[i];
         L = dataptr->slu_Lmat[i];
         U = dataptr->slu_Umat[i];
         perm_c = dataptr->perm_c[i];
         perm_r = dataptr->perm_r[i];
         dCreate_Dense_Matrix(&B, length, 1, rhsbuf, length, DN, _D, GE);
         dCreate_Dense_Matrix(&X, length, 1, solbuf, length, DN, _D, GE);
         dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
                equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
                ferr, berr, &mem_usage, &info);
         for ( j = 0; j < length; j++ ) 
            /* if ( indptr[j] < inlen ) x[indptr[j]] += solbuf[j];*/
if ( indptr[j] < inlen ) x[indptr[j]] = solbuf[j];
         Destroy_SuperMatrix_Store(&B);
         Destroy_SuperMatrix_Store(&X);
      }
   }

   /* --------------------------------------------------------- */
   /* clean up                                                  */
   /* --------------------------------------------------------- */

   if (ntimes > 1) ML_free(xbuffer);
   ML_free( rhsbuf );
   ML_free( solbuf );
   ML_free( dbuffer );
   ML_free( etree );
   SUPERLU_FREE (R);
   SUPERLU_FREE (C);
   SUPERLU_FREE (ferr);
   SUPERLU_FREE (berr);
   return 0;
#else
   printf("ML_Smoother_VBlockAdditiveSchwarz : not available.\n");
   ML_avoid_unused_param( (void *) sm);
   ML_avoid_unused_param( (void *) &inlen);
   ML_avoid_unused_param( (void *) x);
   ML_avoid_unused_param( (void *) &outlen);
   ML_avoid_unused_param( (void *) rhs);
   exit(1);
   return 1;
#endif
}

/* ************************************************************************* */
/* variable block multiplicative Schwarz                                     */
/* ------------------------------------------------------------------------- */

int ML_Smoother_VBlockMultiplicativeSchwarz(ML_Smoother *sm, int inlen, double x[],
                                      int outlen, double rhs[])
{
#ifdef SUPERLU
   int                i, j, k, m, extNrows, nblocks, length, *indptr, ntimes;
   int                index, *blk_size, **blk_indices, max_blk_size;
   int                **aux_bmat_ia, **aux_bmat_ja;
   double             *dbuffer, **aux_bmat_aa, *rhsbuf, *solbuf, *xbuffer=NULL;
   ML_Comm            *comm;
   ML_Operator        *Amat;
   ML_Smoother        *smooth_ptr;
   ML_Sm_Schwarz_Data *dataptr;
   ML_CommInfoOP      *getrow_comm;
   int                info, *perm_r, *perm_c, *etree, panel_size, lwork;
   double             *R, *C, *ferr, *berr, rpg, rcond, dtemp;
   char               fact[1], equed[1], trans[1], refact[1];
   void               *work=NULL;
   SuperMatrix        *A, *L, *U, B, X;
   factor_param_t     iparam;
   mem_usage_t        mem_usage;

   /* --------------------------------------------------------- */
   /* fetch parameters and check                                */
   /* --------------------------------------------------------- */

   smooth_ptr = (ML_Smoother *) sm;
   comm       = smooth_ptr->my_level->comm;
   Amat       = smooth_ptr->my_level->Amat;
   dataptr    = (ML_Sm_Schwarz_Data *) smooth_ptr->smoother->data;
   ntimes     = smooth_ptr->ntimes;

   if (Amat->getrow->ML_id == ML_EMPTY) 
      pr_error("Error(ML_Smoother_MultiplicativeSchwarz): Need getrow()\n");
   if (Amat->getrow->post_comm != NULL)
      pr_error("Post communication not implemented for MultiplicativeSchwarz\n");
   if ( dataptr == NULL ) 
      pr_error("Error(MultiplicativeSchwarz): Need dataptr\n");
   getrow_comm= Amat->getrow->pre_comm;

   extNrows    = dataptr->Nrows;
   nblocks     = dataptr->nblocks;
   blk_indices = dataptr->blk_indices;
   blk_size    = dataptr->blk_size;
   aux_bmat_ia = dataptr->aux_bmat_ia;
   aux_bmat_ja = dataptr->aux_bmat_ja;
   aux_bmat_aa = dataptr->aux_bmat_aa;
   max_blk_size = 0;
   for ( i = 0; i < nblocks; i++ ) 
      max_blk_size = (blk_size[i] > max_blk_size) ? blk_size[i] : max_blk_size;

   /* --------------------------------------------------------- */
   /* communicate the rhs and put into dbuffer                  */
   /* --------------------------------------------------------- */

   dbuffer = (double *) ML_allocate(extNrows * sizeof(double));
   for ( i = 0; i < outlen; i++ ) dbuffer[i] = rhs[i];
   for ( i = outlen; i < extNrows; i++ ) dbuffer[i] = 0.0;
   for ( i = 0; i < inlen;  i++ ) x[i] = 0.0;

   if (extNrows > outlen && getrow_comm != NULL)
      ML_exchange_bdry(dbuffer,getrow_comm,inlen,comm,ML_OVERWRITE,NULL);

   /* --------------------------------------------------------- */
   /* set up for SuperLU solves                                 */
   /* --------------------------------------------------------- */

   rhsbuf = (double *) ML_allocate(max_blk_size * sizeof(double));
   solbuf = (double *) ML_allocate(max_blk_size * sizeof(double));
   panel_size               = sp_ienv(1);
   iparam.panel_size        = panel_size;
   iparam.relax             = sp_ienv(2);
   iparam.diag_pivot_thresh = 1.0;
   iparam.drop_tol          = -1;
   lwork                    = 0;
   *fact                    = 'F';
   *equed                   = 'N';
   *trans                   = 'N';
   *refact                  = 'N';
   R     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   C     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   ferr  = (double *) SUPERLU_MALLOC(sizeof(double));
   berr  = (double *) SUPERLU_MALLOC(sizeof(double));
   etree = (int *) ML_allocate( max_blk_size * sizeof(int) );

   /* --------------------------------------------------------- */
   /* the first pass                                            */
   /* --------------------------------------------------------- */

   for ( i = 0; i < nblocks; i++ )
   {
      indptr = blk_indices[i];
      length = blk_size[i];
      for ( j = 0; j < length; j++ ) rhsbuf[j] = dbuffer[indptr[j]];
      A = dataptr->slu_Amat[i];
      L = dataptr->slu_Lmat[i];
      U = dataptr->slu_Umat[i];
      perm_c = dataptr->perm_c[i];
      perm_r = dataptr->perm_r[i];
      dCreate_Dense_Matrix(&B, length, 1, rhsbuf, length, DN, _D, GE);
      dCreate_Dense_Matrix(&X, length, 1, solbuf, length, DN, _D, GE);
      dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
             equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
             ferr, berr, &mem_usage, &info);
      for ( j = 0; j < length; j++ ) 
         /* if ( indptr[j] < inlen ) x[indptr[j]] += solbuf[j]; */
if ( indptr[j] < inlen ) x[indptr[j]] = solbuf[j];
      Destroy_SuperMatrix_Store(&B);
      Destroy_SuperMatrix_Store(&X);
   }

   if (ntimes > 1) xbuffer = (double *) ML_allocate(extNrows * sizeof(double));

   for ( m = 1; m < ntimes; m++ )
   {
      for ( i = 0; i < inlen; i++ ) xbuffer[i] = x[i];
      if (extNrows > outlen && getrow_comm != NULL)
         ML_exchange_bdry(xbuffer,getrow_comm,inlen,comm,ML_OVERWRITE,NULL);

      for ( i = 0; i < nblocks; i++ )
      {
         indptr = blk_indices[i];
         length = blk_size[i];
         for ( j = 0; j < length; j++ )
         {
            dtemp = dbuffer[indptr[j]];
            for ( k = aux_bmat_ia[i][j]; k < aux_bmat_ia[i][j+1]; k++ )
            {
               index = aux_bmat_ja[i][k];
               if (index < inlen) dtemp -= (aux_bmat_aa[i][k] * x[index]); 
               else               dtemp -= (aux_bmat_aa[i][k] * xbuffer[index]); 
            }
            rhsbuf[j] = dtemp;
         }
         A = dataptr->slu_Amat[i];
         L = dataptr->slu_Lmat[i];
         U = dataptr->slu_Umat[i];
         perm_c = dataptr->perm_c[i];
         perm_r = dataptr->perm_r[i];
         dCreate_Dense_Matrix(&B, length, 1, rhsbuf, length, DN, _D, GE);
         dCreate_Dense_Matrix(&X, length, 1, solbuf, length, DN, _D, GE);
         dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
                equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
                ferr, berr, &mem_usage, &info);
         for ( j = 0; j < length; j++ ) 
            /* if ( indptr[j] < inlen ) x[indptr[j]] += solbuf[j];*/
if ( indptr[j] < inlen ) x[indptr[j]] = solbuf[j];
         Destroy_SuperMatrix_Store(&B);
         Destroy_SuperMatrix_Store(&X);
      }
   }

   /* --------------------------------------------------------- */
   /* clean up                                                  */
   /* --------------------------------------------------------- */

   if (ntimes > 1) ML_free(xbuffer);
   ML_free( rhsbuf );
   ML_free( solbuf );
   ML_free( dbuffer );
   ML_free( etree );
   SUPERLU_FREE (R);
   SUPERLU_FREE (C);
   SUPERLU_FREE (ferr);
   SUPERLU_FREE (berr);
   return 0;
#else
   printf("ML_Smoother_VBlockMultiplicativeSchwarz : not available.\n");
   ML_avoid_unused_param( (void *) sm);
   ML_avoid_unused_param( (void *) &inlen);
   ML_avoid_unused_param( (void *) x);
   ML_avoid_unused_param( (void *) &outlen);
   ML_avoid_unused_param( (void *) rhs);
   exit(1);
   return 1;
#endif
}

/* ************************************************************************* */
/* Constructor for ML_Sm_Schwarz_Data                                        */
/* ************************************************************************* */

int ML_Smoother_Create_Schwarz_Data(ML_Sm_Schwarz_Data **data)
{
   ML_Sm_Schwarz_Data *ml_data;

   ML_memory_alloc((void**) data, sizeof(ML_Sm_Schwarz_Data), "SMI");
   ml_data = (ML_Sm_Schwarz_Data *) (*data);
   ml_data->bmat_ia      = NULL;
   ml_data->bmat_ja      = NULL;
   ml_data->bmat_aa      = NULL;
   ml_data->aux_bmat_ia  = NULL;
   ml_data->aux_bmat_ja  = NULL;
   ml_data->aux_bmat_aa  = NULL;
   ml_data->getrow_comm  = NULL;
   ml_data->Nrows        = 0;
   ml_data->blk_size     = NULL;
   ml_data->blk_info     = NULL;
   ml_data->blk_indices  = NULL;
#ifdef SUPERLU
   ml_data->slu_Amat     = NULL;
   ml_data->slu_Lmat     = NULL;
   ml_data->slu_Umat     = NULL;
#endif 
   ml_data->perm_c       = NULL;
   ml_data->perm_r       = NULL;
   return 0;
}

/* ************************************************************************* */
/* Destructor for ML_Sm_Schwarz_Data                                         */
/* ************************************************************************* */

void ML_Smoother_Destroy_Schwarz_Data(void *data)
{
   int                i;
   ML_Sm_Schwarz_Data *ml_data;
#ifdef SUPERLU
   SuperMatrix        *A, *L, *U;
#endif

   ml_data = (ML_Sm_Schwarz_Data *) data;
   if ( ml_data->bmat_ia  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->bmat_ia[i]);
      ML_free(ml_data->bmat_ia);
   }
   if ( ml_data->bmat_ja  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->bmat_ja[i]);
      ML_free(ml_data->bmat_ja);
   }
   if ( ml_data->bmat_aa  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->bmat_aa[i]);
      ML_free(ml_data->bmat_aa);
   }
   if ( ml_data->aux_bmat_ia  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->aux_bmat_ia[i]);
      ML_free(ml_data->aux_bmat_ia);
   }
   if ( ml_data->aux_bmat_ja  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->aux_bmat_ja[i]);
      ML_free(ml_data->aux_bmat_ja);
   }
   if ( ml_data->aux_bmat_aa  != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ ) ML_free(ml_data->aux_bmat_aa[i]);
      ML_free(ml_data->aux_bmat_aa);
   }
   if ( ml_data->blk_size != NULL ) ML_free(ml_data->blk_size);
   if ( ml_data->blk_info != NULL ) ML_free(ml_data->blk_info);
   if ( ml_data->blk_indices != NULL ) 
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
         if ( ml_data->blk_indices[i] != NULL ) 
            ML_free( ml_data->blk_indices[i] );
   }
#ifdef SUPERLU
   if ( ml_data->slu_Amat != NULL )
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
      {
         A = ml_data->slu_Amat[i];
         if ( A != NULL )
         {
            SUPERLU_FREE( A->Store );
            ML_free(A);
         }
      }
      ML_free( ml_data->slu_Amat );
   }
   if ( ml_data->slu_Lmat != NULL )
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
      {
         L = ml_data->slu_Lmat[i];
         if ( L != NULL )
         {
            Destroy_SuperNode_Matrix(L);
            ML_free(L);
         }
      }
      ML_free( ml_data->slu_Lmat );
   }
   if ( ml_data->slu_Umat != NULL )
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
      {
         U = ml_data->slu_Umat[i];
         if ( U != NULL )
         {
            SUPERLU_FREE( ((NRformat *) U->Store)->colind);
            SUPERLU_FREE( ((NRformat *) U->Store)->rowptr);
            SUPERLU_FREE( ((NRformat *) U->Store)->nzval);
            SUPERLU_FREE( U->Store );
            ML_free(U);
         }
      }
      ML_free( ml_data->slu_Umat );
   }
   if ( ml_data->perm_c != NULL )
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
         if ( ml_data->perm_c[i] ) ML_free(ml_data->perm_c[i]);
      ML_free( ml_data->perm_c );
   }
   if ( ml_data->perm_r != NULL )
   {
      for ( i = 0; i < ml_data->nblocks; i++ )
         if ( ml_data->perm_r[i] ) ML_free(ml_data->perm_r[i]);
      ML_free( ml_data->perm_r );
   }
#endif
   ML_memory_free( (void **) &ml_data);
}

/*****************************************************************************/
/* function for setting up variable block overlapped Schwarz                 */
/*****************************************************************************/

int ML_Smoother_VBlockSchwarzDecomposition(ML_Sm_Schwarz_Data *data, 
             ML_Operator *Amat, ML_Comm *comm, int total_recv_leng, 
             int *recv_lengths, int *ext_ja, double *ext_aa, int *map, 
             int *map2, int Noffset)
{
#ifdef SUPERLU
   int                i, j, k, **bmat_ia, **bmat_ja, allocated_space;
   int                *blk_size, index, **blk_indices, **aux_bmat_ia;
   int                offset, nnz, Nrows, extNrows, **aux_bmat_ja;
   int                mypid, *tmp_blk_leng, *cols, blknum, ncnt, *blkinfo;
   int                rownum, rowleng, nblocks, col_ind, init_size, aux_nnz;
   int                *tmp_indices, cur_off_row, max_blk_size, *mat_ia, *mat_ja;
   double             **bmat_aa, *vals, *mat_aa, **aux_bmat_aa;
   ML_Sm_Schwarz_Data *schwarz_ptr;
   int                info, *perm_r, *perm_c, permc_spec, *etree, panel_size;
   int                lwork, nrows;
   double             *R, *C, *ferr, *berr, rpg, rcond, *trhs, *tsol;
   char               fact[1], equed[1], trans[1], refact[1];
   void               *work=NULL;
   SuperMatrix        *A, *L, *U, B, X;
   factor_param_t     iparam;
   mem_usage_t        mem_usage;

   /* ---------------------------------------------------------- */
   /* fetch Schwarz parameters                                   */
   /* ---------------------------------------------------------- */

   mypid       = comm->ML_mypid;
   schwarz_ptr = (ML_Sm_Schwarz_Data *) data;
   Nrows       = Amat->outvec_leng;
   extNrows    = Nrows + total_recv_leng;
   schwarz_ptr->Nrows = extNrows;
   blkinfo     = schwarz_ptr->blk_info;
   nblocks     = schwarz_ptr->nblocks;

   /* ---------------------------------------------------------- */
   /* adjust the off-processor row data                          */
   /* ---------------------------------------------------------- */

   offset = 0;
   for ( i = 0; i < total_recv_leng; i++ )
   {
      for ( j = offset; j < offset+recv_lengths[i]; j++ ) 
      {
         index = ext_ja[j];
         if ( index >= Noffset && index < Noffset+Nrows )
            ext_ja[j] = index - Noffset; 
         else
         {
            col_ind = ML_sorted_search(index, extNrows-Nrows, map); 
            if ( col_ind >= 0 ) ext_ja[j] = map2[col_ind] + Nrows;
            else                ext_ja[j] = -1;
         }
      }
      offset += recv_lengths[i];
   }

   /* ---------------------------------------------------------- */
   /* compose the initial blk_size information                   */
   /* ---------------------------------------------------------- */

   schwarz_ptr->blk_size = (int *) ML_allocate(nblocks * sizeof(int) );
   schwarz_ptr->blk_indices = (int **) ML_allocate(nblocks * sizeof(int*) );
   blk_indices  = schwarz_ptr->blk_indices;
   blk_size     = schwarz_ptr->blk_size;
   tmp_blk_leng = (int *) ML_allocate(nblocks * sizeof(int) );
   for ( i = 0; i < nblocks; i++ ) blk_size[i] = 0; 
   for ( i = 0; i < Nrows; i++ )   blk_size[blkinfo[i]]++;
   for ( i = 0; i < nblocks; i++ ) 
   {
      if ( blk_size[i] == 0 )
      {
         printf("%4d : SchwarzDecomposition - block %d is empty\n",mypid,i);
         exit(1);
      }
   }
   for ( i = 0; i < nblocks; i++ ) 
   {
      tmp_blk_leng[i] = blk_size[i] * blk_size[i] + 5;
      blk_indices[i] = (int *) ML_allocate(tmp_blk_leng[i] * sizeof(int));
   }
   for ( i = 0; i < nblocks; i++ ) blk_size[i] = 0;
   for ( i = 0; i < Nrows; i++ ) 
   {
      blknum = blkinfo[i];
      index = blk_size[blknum]++;
      blk_indices[blknum][index] = i;
   }

   /* ---------------------------------------------------------- */
   /* now extend the each block for the overlap                  */
   /* (at the end blk_indices and bli_size contains the info)    */
   /* ---------------------------------------------------------- */

   allocated_space = extNrows;
   vals = (double *) ML_allocate(allocated_space * sizeof(double));
   cols = (int *)    ML_allocate(allocated_space * sizeof(int));
   max_blk_size = 0;
   for ( i = 0; i < nblocks; i++ ) 
   {
      init_size = blk_size[i];
      for ( j = 0; j < init_size; j++ ) 
      {
         rownum = blk_indices[i][j];
         ML_get_matrix_row(Amat,1,&rownum,&allocated_space,&cols,&vals,&rowleng,0);
         if ( blk_size[i] + rowleng > tmp_blk_leng[i] )
         {
            tmp_indices = blk_indices[i];
            tmp_blk_leng[i] = 2 * ( blk_size[i] + rowleng ) + 2;
            blk_indices[i] = (int *) ML_allocate(tmp_blk_leng[i] * sizeof(int));
            for (k = 0; k < blk_size[i]; k++) blk_indices[i][k] = tmp_indices[k]; 
            ML_free( tmp_indices );
         }   
         for ( k = 0; k < rowleng; k++ ) 
         {
            col_ind = cols[k];
            blk_indices[i][blk_size[i]++] = col_ind;
         }
      } 
      ML_az_sort(blk_indices[i], blk_size[i], NULL, NULL);
      ncnt = 0;
      for ( j = 1; j < blk_size[i]; j++ ) 
         if ( blk_indices[i][j] != blk_indices[i][ncnt] )
           blk_indices[i][++ncnt] = blk_indices[i][j];
      blk_size[i] = ncnt + 1;
      if ( blk_size[i] > max_blk_size ) max_blk_size = blk_size[i];
   }

   /* ---------------------------------------------------------- */
   /* compute the memory requirements for each block             */
   /* ---------------------------------------------------------- */

   schwarz_ptr->bmat_ia = (int **)    ML_allocate(nblocks * sizeof(int*) );
   schwarz_ptr->bmat_ja = (int **)    ML_allocate(nblocks * sizeof(int*) );
   schwarz_ptr->bmat_aa = (double **) ML_allocate(nblocks * sizeof(double*) );
   bmat_ia = schwarz_ptr->bmat_ia;
   bmat_ja = schwarz_ptr->bmat_ja;
   bmat_aa = schwarz_ptr->bmat_aa;
   schwarz_ptr->aux_bmat_ia = (int **)    ML_allocate(nblocks * sizeof(int*) );
   schwarz_ptr->aux_bmat_ja = (int **)    ML_allocate(nblocks * sizeof(int*) );
   schwarz_ptr->aux_bmat_aa = (double **) ML_allocate(nblocks * sizeof(double*) );
   aux_bmat_ia = schwarz_ptr->aux_bmat_ia;
   aux_bmat_ja = schwarz_ptr->aux_bmat_ja;
   aux_bmat_aa = schwarz_ptr->aux_bmat_aa;

   for ( i = 0; i < nblocks; i++ ) 
   {
      nnz = aux_nnz = offset = cur_off_row = 0;
      for ( j = 0; j < blk_size[i]; j++ ) 
      {
         rownum = blk_indices[i][j];
         if ( rownum < Nrows )
            ML_get_matrix_row(Amat,1,&rownum,&allocated_space,&cols,
                              &vals,&rowleng,0);
         else 
         {
            for ( k = cur_off_row; k < rownum-Nrows; k++ ) 
               offset += recv_lengths[k]; 
            cur_off_row = rownum - Nrows;
            rowleng = 0;
            for ( k = offset; k < offset+recv_lengths[cur_off_row]; k++ ) 
               if ( ext_ja[k] != -1 ) cols[rowleng++] = ext_ja[k];
         }
         for ( k = 0; k < rowleng; k++ )
         {
            index = ML_find_index( cols[k], blk_indices[i], blk_size[i]);
            if ( index >= 0 ) nnz++;
            else              aux_nnz++;
         }
      }
      bmat_ia[i] = (int *)    ML_allocate( (blk_size[i] + 1) * sizeof(int));
      bmat_ja[i] = (int *)    ML_allocate( nnz * sizeof(int));
      bmat_aa[i] = (double *) ML_allocate( nnz * sizeof(double));
      aux_bmat_ia[i] = (int *)    ML_allocate( (blk_size[i] + 1) * sizeof(int));
      aux_bmat_ja[i] = (int *)    ML_allocate( aux_nnz * sizeof(int));
      aux_bmat_aa[i] = (double *) ML_allocate( aux_nnz * sizeof(double));
   }

   /* ---------------------------------------------------------- */
   /* load the submatrices                                       */
   /* ---------------------------------------------------------- */

   for ( i = 0; i < nblocks; i++ ) 
   {
      nnz = aux_nnz     = offset = cur_off_row = 0;
      bmat_ia[i][0]     = 0;
      aux_bmat_ia[i][0] = 0;

      for ( j = 0; j < blk_size[i]; j++ ) 
      {
         rownum = blk_indices[i][j];
         if ( rownum < Nrows )
            ML_get_matrix_row(Amat,1,&rownum,&allocated_space,&cols,
                              &vals,&rowleng,0);
         else 
         {
            for ( k = cur_off_row; k < rownum-Nrows; k++ ) 
            {
               offset += recv_lengths[k]; 
            }
            cur_off_row = rownum - Nrows;
            rowleng = 0;
            for ( k = offset; k < offset+recv_lengths[cur_off_row]; k++ ) 
            {
               if ( ext_ja[k] != -1 ) 
               {
                  cols[rowleng] = ext_ja[k];
                  vals[rowleng++] = ext_aa[k];
               }
            }
         }
         for ( k = 0; k < rowleng; k++ )
         {
            index = ML_find_index( cols[k], blk_indices[i], blk_size[i]);
            if ( index >= 0 )
            {
               bmat_ja[i][nnz] = index;
               bmat_aa[i][nnz++] = vals[k];
            }
            else
            {
               aux_bmat_ja[i][aux_nnz] = cols[k];
               aux_bmat_aa[i][aux_nnz++] = vals[k];
            }
         }
         bmat_ia[i][j+1] = nnz;
         aux_bmat_ia[i][j+1] = aux_nnz;
      } 
   } 

   /* ---------------------------------------------------------- */
   /* call SuperLU to perform decomposition                      */
   /* ---------------------------------------------------------- */

   schwarz_ptr->slu_Amat = (SuperMatrix **) ML_allocate(nblocks*sizeof(SuperMatrix*));
   schwarz_ptr->slu_Lmat = (SuperMatrix **) ML_allocate(nblocks*sizeof(SuperMatrix*));
   schwarz_ptr->slu_Umat = (SuperMatrix **) ML_allocate(nblocks*sizeof(SuperMatrix*));
   schwarz_ptr->perm_r   = (int **) ML_allocate(nblocks*sizeof(int*));
   schwarz_ptr->perm_c   = (int **) ML_allocate(nblocks*sizeof(int*));
   etree = (int *) ML_allocate( max_blk_size * sizeof(int) );
   R     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   C     = (double *) SUPERLU_MALLOC(max_blk_size * sizeof(double));
   ferr  = (double *) SUPERLU_MALLOC(sizeof(double));
   berr  = (double *) SUPERLU_MALLOC(sizeof(double));
   tsol  = (double *) ML_allocate( max_blk_size * sizeof(double) );
   trhs  = (double *) ML_allocate( max_blk_size * sizeof(double) );
   for ( i = 0; i < max_blk_size; i ++ ) trhs[i] = 0.0;

   for ( i = 0; i < nblocks; i ++ )
   {
      schwarz_ptr->slu_Amat[i] = (SuperMatrix *) ML_allocate(sizeof(SuperMatrix));
      schwarz_ptr->slu_Lmat[i] = (SuperMatrix *) ML_allocate(sizeof(SuperMatrix));
      schwarz_ptr->slu_Umat[i] = (SuperMatrix *) ML_allocate(sizeof(SuperMatrix));
      A = schwarz_ptr->slu_Amat[i];
      L = schwarz_ptr->slu_Lmat[i];
      U = schwarz_ptr->slu_Umat[i];
      nrows  = blk_size[i];
      mat_ia = schwarz_ptr->bmat_ia[i];
      mat_ja = schwarz_ptr->bmat_ja[i];
      mat_aa = schwarz_ptr->bmat_aa[i];
      nnz    = mat_ia[nrows];
      dCreate_CompRow_Matrix(A,nrows,nrows,nnz,mat_aa,mat_ja,mat_ia,NR,_D,GE);
      schwarz_ptr->perm_r[i] = (int *) ML_allocate(nrows*sizeof(int));
      schwarz_ptr->perm_c[i] = (int *) ML_allocate(2*nrows*sizeof(int));
      perm_r = schwarz_ptr->perm_r[i];
      perm_c = schwarz_ptr->perm_c[i];
      permc_spec = 0;
      get_perm_c(permc_spec, A, perm_c);
      panel_size               = sp_ienv(1);
      iparam.panel_size        = panel_size;
      iparam.relax             = sp_ienv(2);
      iparam.diag_pivot_thresh = 1.0;
      iparam.drop_tol          = -1;
      lwork                    = 0;
      *fact                    = 'N';
      *equed                   = 'N';
      *trans                   = 'N';
      *refact                  = 'N';
      dCreate_Dense_Matrix(&B, nrows, 1, trhs, nrows, DN, _D, GE);
      dCreate_Dense_Matrix(&X, nrows, 1, tsol, nrows, DN, _D, GE);

      dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
             equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
             ferr, berr, &mem_usage, &info);
#ifdef ML_DEBUG_SMOOTHER
printf("Block %6d(%6d) : cond no. = %e\n", i, nblocks, 1.0 / rcond);
if ( info != 0 && info != (nrows+1) )
{
   for ( j = 0; j < nrows; j++ )
   {
      for ( k = mat_ia[j]; k < mat_ia[j+1]; k++ )
         printf("Block %4d : in data %6d(%6d) %6d = %e\n",i,j,
                nrows,blk_indices[i][mat_ja[k]],mat_aa[k]);
      for ( k = aux_bmat_ia[i][j]; k < aux_bmat_ia[i][j+1]; k++ )
         printf("Block %4d : offdata %6d(%6d) %6d = %e\n",i,j,nrows,
                aux_bmat_ja[i][k],aux_bmat_aa[i][k]);
   }
}
#endif
      Destroy_SuperMatrix_Store(&B);
      Destroy_SuperMatrix_Store(&X);
   }

   /* ---------------------------------------------------------- */
   /* clean up                                                   */
   /* ---------------------------------------------------------- */

   SUPERLU_FREE (R);
   SUPERLU_FREE (C);
   SUPERLU_FREE (ferr);
   SUPERLU_FREE (berr);
   ML_free(etree);
   ML_free(trhs);
   ML_free(tsol);
   ML_free(vals);
   ML_free(cols);
   ML_free(tmp_blk_leng);
   return 0;
#else
   printf("ML_Smoother_VBlockSchwarzDecomposition : not available.\n");
   ML_avoid_unused_param( (void *) data);
   ML_avoid_unused_param( (void *) Amat);
   ML_avoid_unused_param( (void *) comm);
   ML_avoid_unused_param( (void *) &total_recv_leng);
   ML_avoid_unused_param( (void *) recv_lengths);
   ML_avoid_unused_param( (void *) ext_ja);
   ML_avoid_unused_param( (void *) ext_aa);
   ML_avoid_unused_param( (void *) map);
   ML_avoid_unused_param( (void *) map2);
   ML_avoid_unused_param( (void *) &Noffset);
   exit(1);
   return 1;
#endif
}

/* ------------------------------------------------------------------------- */
/* generate the variable block additive Schwarz smoother                     */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockAdditiveSchwarz(ML *ml , int nl, int pre_or_post,
                                          int ntimes, int length, int *blkinfo)
{
   int                (*fun)(ML_Smoother *, int, double *, int, double *);
   int                total_recv_leng, *recv_lengths, *int_buf, *map, *map2; 
   int                i, maxblk, offset;
   double             *dble_buf;
   ML_Sm_Schwarz_Data *data;
   ML_Operator        *Amat;
   ML_Comm            *comm;
   char               str[80];
	
   /* ---------------------------------------------------------------------- */
   /* check for valid incoming data                                          */
   /* ---------------------------------------------------------------------- */
   if (nl == ML_ALL_LEVELS) { 
      printf("ML_Gen_Smoother_VBlockAdditiveSchwarz: ML_ALL_LEVELS not allowed\n");
      return 1;
   }
   if (nl < 0) {
      printf("ML_Gen_Smoother_VBlockAdditiveSchwarz: cannot set smoother on level %d\n",nl);
      return 1;
   }

   Amat = &(ml->Amat[nl]);
   if ( length != 0 && length != Amat->outvec_leng )
   {
      printf("ML_Gen_Smoother_VBlockAdditiveSchwarz ERROR : invalid length.\n");
      exit(1);
   }

   /* ---------------------------------------------------------------------- */
   /* set the nblock and blk_info data                                       */
   /* ---------------------------------------------------------------------- */

   fun = ML_Smoother_VBlockAdditiveSchwarz;
	
   comm = ml->comm;
   ML_Smoother_Create_Schwarz_Data( &data );
   data->Nrows   = Amat->outvec_leng;
   data->blk_info = (int *) ML_allocate(data->Nrows * sizeof(int));
   if ( blkinfo != NULL && length != 0 )
   {
      for ( i = 0; i < length; i++ ) data->blk_info[i] = blkinfo[i];
      maxblk = 0;
      for ( i = 0; i < length; i++ ) 
         if ( blkinfo[i] > maxblk ) maxblk = blkinfo[i];
      data->nblocks = maxblk + 1;
   }
   else 
   {
      for ( i = 0; i < data->Nrows; i++ ) data->blk_info[i] = i;
      data->nblocks = data->Nrows;
   }

   /* ---------------------------------------------------------------- */
   /* send the lengths of each row to remote processor at the end,     */
   /* additional row information should be given in total_recv_leng,   */
   /* recv_lengths, int_buf, dble_buf                                  */
   /* ---------------------------------------------------------------- */

   ML_Smoother_ComposeOverlappedMatrix(Amat, comm, &total_recv_leng, 
              &recv_lengths, &int_buf, &dble_buf, &map, &map2, &offset);

   /* ---------------------------------------------------------------- */
   /* use the local matrix row and the off-processor rows to compose   */
   /* Schwarz preconditioner                                           */
   /* ---------------------------------------------------------------- */

   ML_Smoother_VBlockSchwarzDecomposition(data,Amat,comm,total_recv_leng,
              recv_lengths, int_buf, dble_buf, map, map2, offset);
   if ( map  != NULL ) ML_free(map);
   if ( map2 != NULL ) ML_free(map2);
   if ( int_buf != NULL ) ML_free(int_buf);
   if ( dble_buf != NULL ) ML_free(dble_buf);
   if ( recv_lengths != NULL ) ML_free(recv_lengths);

   /* ---------------------------------------------------------------- */
   /* set it up as smoother                                            */
   /* ---------------------------------------------------------------- */

   if (pre_or_post == ML_PRESMOOTHER) {
      ml->pre_smoother[nl].data_destroy = ML_Smoother_Destroy_Schwarz_Data;
      sprintf(str,"VBASz_pre%d",nl);
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                        (void *) data, fun, NULL, ntimes, 0.0, str));
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_Schwarz_Data;
      sprintf(str,"VBASz_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, 0.0, str));
   }
   else if (pre_or_post == ML_BOTH) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_Schwarz_Data;
      sprintf(str,"VBASz_pre%d",nl);
      ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                      (void *) data, fun, NULL, ntimes, 0.0, str);
      sprintf(str,"VBASz_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, 0.0, str));
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the variable block additive Schwarz smoother                     */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockMultiplicativeSchwarz(ML *ml , int nl, int pre_or_post,
                                         int ntimes, int length, int *blkinfo )
{
   int                (*fun)(ML_Smoother *, int, double *, int, double *);
   int                total_recv_leng, *recv_lengths, *int_buf, *map, *map2; 
   int                i, maxblk, offset;
   double             *dble_buf;
   ML_Sm_Schwarz_Data *data;
   ML_Operator        *Amat;
   ML_Comm            *comm;
   char               str[80];
	
   /* ---------------------------------------------------------------------- */
   /* check for valid incoming data                                          */
   /* ---------------------------------------------------------------------- */

   if (nl == ML_ALL_LEVELS) { 
      printf("ML_Gen_Smoother_VBlockMultiplicativeSchwarz: ML_ALL_LEVELS not allowed\n");
      return 1;
   }
   if (nl < 0) {
      printf("ML_Gen_Smoother_VBlockMultiplicativeSchwarz: cannot set smoother on level %d\n",nl);
      return 1;
   }
   Amat = &(ml->Amat[nl]);
   if ( length != 0 && length != Amat->outvec_leng )
   {
      printf("ML_Gen_Smoother_VBlockMultiplicativeSchwarz : invalid length.\n");
      exit(1);
   }

   /* ---------------------------------------------------------------------- */
   /* set the nblock and blk_info data                                       */
   /* ---------------------------------------------------------------------- */

   fun = ML_Smoother_VBlockMultiplicativeSchwarz;
	
   comm = ml->comm;
   ML_Smoother_Create_Schwarz_Data( &data );
   data->Nrows   = Amat->outvec_leng;
   data->blk_info = (int *) ML_allocate(data->Nrows * sizeof(int));
   if ( blkinfo != NULL && length != 0 )
   {
      for ( i = 0; i < length; i++ ) data->blk_info[i] = blkinfo[i];
      maxblk = 0;
      for ( i = 0; i < length; i++ ) 
         if ( blkinfo[i] > maxblk ) maxblk = blkinfo[i];
      data->nblocks = maxblk + 1;
   }
   else 
   {
      for ( i = 0; i < data->Nrows; i++ ) data->blk_info[i] = i;
      data->nblocks = data->Nrows;
   }

   /* ---------------------------------------------------------------- */
   /* send the lengths of each row to remote processor at the end,     */
   /* additional row information should be given in total_recv_leng,   */
   /* recv_lengths, int_buf, dble_buf                                  */
   /* ---------------------------------------------------------------- */

   ML_Smoother_ComposeOverlappedMatrix(Amat, comm, &total_recv_leng, 
              &recv_lengths, &int_buf, &dble_buf, &map, &map2, &offset);

   /* ---------------------------------------------------------------- */
   /* use the local matrix row and the off-processor rows to compose   */
   /* Schwarz preconditioner                                           */
   /* ---------------------------------------------------------------- */

   ML_Smoother_VBlockSchwarzDecomposition(data,Amat,comm,total_recv_leng,
              recv_lengths, int_buf, dble_buf, map, map2, offset);
   if ( map  != NULL ) ML_free(map);
   if ( map2 != NULL ) ML_free(map2);
   if ( int_buf != NULL ) ML_free(int_buf);
   if ( dble_buf != NULL ) ML_free(dble_buf);
   if ( recv_lengths != NULL ) ML_free(recv_lengths);

   /* ---------------------------------------------------------------- */
   /* set it up as smoother                                            */
   /* ---------------------------------------------------------------- */

   if (pre_or_post == ML_PRESMOOTHER) {
      ml->pre_smoother[nl].data_destroy = ML_Smoother_Destroy_Schwarz_Data;
      sprintf(str,"VBMSz_pre%d",nl);
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                        (void *) data, fun, NULL, ntimes, 0.0, str));
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_Schwarz_Data;
      sprintf(str,"VBMSz_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, 0.0, str));
   }
   else if (pre_or_post == ML_BOTH) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_Schwarz_Data;
      sprintf(str,"VBMSz_pre%d",nl);
      ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                        (void *) data, fun, NULL, ntimes, 0.0, str);
      sprintf(str,"VBMSz_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, 0.0, str));
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
}





#ifdef WKC
/* WKC  The double * is actually an Epetra_MultiVector * */

/* NO NEED TO BLOCK, SOLVING A #proc x #proc system at this point, 
   communication costs vastly dominates the calculation.  */

int ML_SuperLU_Solve_WKC(ML_Solver *vsolver,int ilen,double *pep_x,int olen,double *pep_rhs)
{



#ifdef SUPERLU
double **pp_x , **pp_rhs;
Epetra_MultiVector &ep_x ( *(Epetra_MultiVector *)pep_x );
Epetra_MultiVector &ep_rhs ( *(Epetra_MultiVector *)pep_rhs );
ep_x.ExtractView ( &pp_x );
ep_rhs.ExtractView ( &pp_rhs );
for ( int KK = 0 ; KK != ep_x.NumVectors() ; KK++ ) {
   double *x = pp_x[KK];
   double *rhs = pp_rhs[KK];


   int            i, n, info, flag, *perm_r, *perm_c, permc_spec;
   int            N_local, offset, *etree, panel_size, lwork;
   double         *local_x, *local_rhs, *R, *C, *ferr, *berr;
   double         rpg, rcond;
   char           fact[1], equed[1], trans[1], refact[1];
   void           *work=NULL;
   ML_Comm        *comm;
   ML_Solver      *solver;
   SuperMatrix    *A, *L, *U, B, X;
   factor_param_t iparam;
   mem_usage_t    mem_usage;

#ifdef ML_TIMING_DETAILED
   double         t0;
#endif

   /* ------------------------------------------------------------- */
   /* fetch the sparse matrix and other parameters                  */
   /* ------------------------------------------------------------- */

   if ( ilen != olen ) 
   {
      printf("ML_SuperLU_Solve error : lengths not matched.\n");
      exit(1);
   }
   solver   = (ML_Solver *) vsolver;

   comm     = (ML_Comm *) solver->void_params1;
   N_local  = (int) solver->dble_params1[0];
   offset   = (int) solver->dble_params1[1];
   n        = (int) solver->dble_params1[2];
   flag     = solver->reuse_flag;

   /* ------------------------------------------------------------- */
   /* if factorization has not been done, allocate space for it     */
   /* ------------------------------------------------------------- */

   if ( flag == -999 ) 
   {
      A = (SuperMatrix *) solver->Mat1;
      if ( A != NULL )
      {
         ML_memory_free((void**)&(((NRformat *) A->Store)->colind));
         ML_memory_free((void**)&(((NRformat *) A->Store)->rowptr));
         ML_memory_free((void**)&(((NRformat *) A->Store)->nzval));
/*WKC         ML_memory_free((void*)&(((NRformat *) A->Store)->colind));
  WKC         ML_memory_free((void*)&(((NRformat *) A->Store)->rowptr));
  WKC         ML_memory_free((void*)&(((NRformat *) A->Store)->nzval)); */
         SUPERLU_FREE( ((NRformat *) A->Store)->nzval);
         SUPERLU_FREE( A->Store );
         ML_memory_free((void**) &A);
      }
      solver->Mat1 = NULL;
      L = (SuperMatrix *) solver->Mat2;
      if ( L != NULL ) 
      {
         Destroy_SuperNode_Matrix(L);
         ML_memory_free((void**) &L);
      }
      solver->Mat2 = NULL;
      U = (SuperMatrix *) solver->Mat3;
      if ( U != NULL )
      {
         SUPERLU_FREE( ((NRformat *) U->Store)->colind);
         SUPERLU_FREE( ((NRformat *) U->Store)->rowptr);
         SUPERLU_FREE( ((NRformat *) U->Store)->nzval);
         SUPERLU_FREE( U->Store );
         ML_memory_free((void**) &U);
      }
      solver->Mat3 = NULL;
      perm_r = (int *) solver->int_params1;
      if (perm_r != NULL) ML_memory_free((void**) &(solver->int_params1));
      solver->int_params1 = NULL;
      perm_c = (int *) solver->int_params2;
      if (perm_c != NULL) ML_memory_free((void**) &(solver->int_params2));
      solver->int_params2 = NULL;
#ifdef ML_TIMING_DETAILED
      if ( comm->ML_mypid == 0 )
         printf("Total SuperLU solve time = %e\n", solver->dble_data);
#endif
      return 0;
   } 
   else if ( flag == 0 ) 
   {
      A = (SuperMatrix *) solver->Mat1;
      ML_memory_alloc((void**) &perm_c, n * sizeof(int), "LU2" );
      ML_memory_alloc((void**) &perm_r, 2 * n * sizeof(int), "LU3" );
      solver->int_params1 = perm_r;
      solver->int_params2 = perm_c;
      permc_spec = 0;
      get_perm_c(permc_spec, A, perm_c);
      ML_memory_alloc((void**) &L, sizeof(SuperMatrix), "LU4" );
      ML_memory_alloc((void**) &U, sizeof(SuperMatrix), "LU5" );
      solver->Mat2 = (void *) L;
      solver->Mat3 = (void *) U;
   } 
   else 
   {
      A = (SuperMatrix *) solver->Mat1;
      L = (SuperMatrix *) solver->Mat2;
      U = (SuperMatrix *) solver->Mat3;
      perm_r = (int *) solver->int_params1;
      perm_c = (int *) solver->int_params2;
   }

   /* ------------------------------------------------------------- */
   /* gather from all processors the complete right hand side       */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &local_rhs, n*sizeof(double),"LU1" );
   ML_memory_alloc((void**) &local_x,   n*sizeof(double),"LU2" );
   for ( i = 0; i < N_local; i++ ) local_rhs[i] = rhs[i];
   i = N_local;
   ML_Comm_GappendDouble((ML_Comm *) comm, local_rhs, &i, n);

   /* ------------------------------------------------------------- */
   /* create right hand side vector conforming to SuperLU format    */
   /* ------------------------------------------------------------- */

   dCreate_Dense_Matrix(&B, n, 1, local_rhs, n, DN, _D, GE);
   dCreate_Dense_Matrix(&X, n, 1, local_x,   n, DN, _D, GE);

   /* ------------------------------------------------------------- */
   /* perform LU decomposition and then solve                       */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &etree,   n*sizeof(int),"LU3" );
   panel_size               = sp_ienv(1);
   iparam.panel_size        = panel_size;
   iparam.relax             = sp_ienv(2);
   iparam.diag_pivot_thresh = 1.0;
   iparam.drop_tol          = -1;
   lwork                    = 0;
   if ( flag == 0 )   *fact = 'N';
   else               *fact = 'F';
   *equed                   = 'N';
   *trans                   = 'N';
   *refact                  = 'N';
   R    = (double *) SUPERLU_MALLOC(n * sizeof(double));
   C    = (double *) SUPERLU_MALLOC(n * sizeof(double));
   ferr = (double *) SUPERLU_MALLOC(sizeof(double));
   berr = (double *) SUPERLU_MALLOC(sizeof(double));

#ifdef ML_TIMING_DETAILED
   t0 = GetClock();
#endif

   dgssvx(fact, trans, refact, A, &iparam, perm_c, perm_r, etree,
          equed, R, C, L, U, work, lwork, &B, &X, &rpg, &rcond,
          ferr, berr, &mem_usage, &info);

#ifdef ML_TIMING_DETAILED
   t0 = GetClock() - t0;
   solver->dble_data += t0;
#endif

   if ( info != 0 && info != n+1 )
   {
      printf("ML_SuperLU_Solve : error coming from dgssvx %d\n", info);
      exit(1);
   } 
   else if ( solver->reuse_flag == 0 )
   {
#ifdef ML_DEBUG_SUPERLU
      if ( rcond != 0.0 && offset == 0 )
         printf("ML_SuperLU_Solve : condition number = %e\n", 1.0/rcond);
      else if ( offset == 0 )
         printf("ML_SuperLU_Solve : Recip. condition number = %e\n", rcond);
#endif
   }

   /* ------------------------------------------------------------- */
   /* extract the local solution sub-vector and then clean up       */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < N_local; i++ ) x[i] = local_x[i+offset];

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   SUPERLU_FREE (R);
   SUPERLU_FREE (C);
   SUPERLU_FREE (ferr);
   SUPERLU_FREE (berr);
   ML_memory_free( (void **) &local_x );
   ML_memory_free( (void **) &local_rhs );
   ML_memory_free( (void **) &etree );
   Destroy_SuperMatrix_Store(&B);
   Destroy_SuperMatrix_Store(&X);
   solver->reuse_flag = 1;
}
#else
#ifdef DSUPERLU
   int                flag, N_local, offset;
   double            *local_rhs;
   ML_Comm           *mlcomm;
   ML_Solver         *solver;
   superlu_options_t  options;
   superlu_options_t *optionsptr;
   SuperMatrix       *A;
   ScalePermstruct_t *ScalePermstruct;
   SuperLUStat_t      stat;
   LUstruct_t        *LUstruct;
   gridinfo_t        *mygrid;
   ML_Lugrid         *lugrid_tiles;
   double             berr;
   int_t              i, n, nprow, npcol, nprocs; /* short integers */
   int                iam, info, nrhs, color, key;
   int                q, g, l, k, mygroup;
   int                stile, mtile, ltile, ntile, tsz, tsz2;
   int_t             *usermap;
   /* In-lining superlu_gridmap */
   MPI_Group mpi_base_group, superlu_grp;
   int                mycol, myrow, j;
   /* heap info arguments */
   int fragments, total_free, largest_free, total_used;

   /* ------------------------------------------------------------- */
   /* fetch the sparse matrix and other parameters                  */
   /* ------------------------------------------------------------- */

   if ( ilen != olen ) 
   {
      printf("ML_SuperLU_Solve error : lengths not matched.\n");
      exit(1);
   }
   solver          = (ML_Solver *) vsolver;
   A               = (SuperMatrix *) solver->Mat1;
   mlcomm          = (ML_Comm *) solver->void_params1;
   iam             = mlcomm->ML_mypid;
   nprocs          = mlcomm->ML_nprocs;
   N_local         = (int) solver->dble_params1[0];
   offset          = (int) solver->dble_params1[1];
   i               = (int) solver->dble_params1[2];
   n               = i;   /* n is a short integer */
   flag            = solver->reuse_flag;
   optionsptr      = &options;
   LUstruct        = (LUstruct_t *) solver->LUspl;
   ScalePermstruct = (ScalePermstruct_t *) solver->PERMspl;

   /* ------------------------------------------------------------- */
   /* if factorization has not been done, allocate space for it     */
   /* Fetching the factorization (flag=1) is unneccessary           */
   /* ------------------------------------------------------------- */

   if ( flag == -999 ) 
   {
     if( iam == 0 )printf("ml_superlu: clean up\n"); /* dmd */

     /* deallocate storage and clean up */
     info = flag;
     if ( A != NULL ) 
     {
        ML_memory_free((void*)&(((NCformat *) A->Store)->rowind));
        ML_memory_free((void*)&(((NCformat *) A->Store)->colptr));
        ML_memory_free((void*)&(((NCformat *) A->Store)->nzval));
        SUPERLU_FREE( A->Store );
        ML_memory_free((void**) &A);
     }
     solver->Mat1 = NULL;
     mygroup = (int ) solver->ML_subgroup;
     lugrid_tiles    = solver->gridtiles;
     mygrid = &((lugrid_tiles[mygroup]).grid);
     Destroy_LU(n, mygrid, LUstruct);
     ScalePermstructFree(ScalePermstruct);
     LUstructFree(LUstruct);
     ML_free(LUstruct);
     ML_free(ScalePermstruct);
     solver->PERMspl = NULL;
     solver->LUspl = NULL;
     superlu_gridexit(mygrid);
     free (lugrid_tiles);
     solver->gridtiles = NULL;
     return 0;
   } 
   else if ( flag == 0 ) 
   {
      ML_SuperLU_Set_Tile(nprocs, &tsz ,&stile, &mtile, &ltile);
      ntile = stile + mtile + ltile;
      tsz2 = tsz * tsz;
      usermap = (int_t *) ML_allocate( tsz2*sizeof(int_t) );
      lugrid_tiles = (ML_Lugrid *) ML_allocate( ntile*sizeof(ML_Lugrid) );
      k = 0;
      for( g=0 ; g < ntile; g++)
      {
         if( g < stile )
         {
            tsz2 = (tsz-1)*(tsz-1);
            nprow = tsz-1;
            npcol = tsz-1;
         }
         else if( g < stile+mtile)
         {
            tsz2 = tsz*(tsz-1);
            nprow = tsz;
            npcol = tsz-1;
         }
         else
         {
            tsz2 = tsz*tsz;
            nprow = tsz;
            npcol = tsz;
         }
         for( l=0; l<tsz2; l++)
         {
            usermap[l] = l+k;
            if( iam == l+k ) mygroup = g;
         }
         k = k + tsz2;
         /* in-lining
          * superlu_gridmap( MPI_COMM_WORLD, 
          * nprow, npcol, usermap, nprow, &((lugrid_tiles[g]).grid));
          */
         (lugrid_tiles[g]).grid.nprow = nprow;
         (lugrid_tiles[g]).grid.npcol = npcol;
         MPI_Comm_group( MPI_COMM_WORLD, &mpi_base_group );
         MPI_Group_incl( mpi_base_group, tsz2, usermap, &superlu_grp );
         MPI_Comm_create(MPI_COMM_WORLD,superlu_grp,&(lugrid_tiles[g].grid.comm));
         if ( lugrid_tiles[g].grid.comm == MPI_COMM_NULL ) 
         {
            lugrid_tiles[g].grid.comm = MPI_COMM_WORLD;
            lugrid_tiles[g].grid.iam  = iam;
         } 
         else 
         {   /* mygroup=g and iam%mtile=npcol*myrow + mycol */
            MPI_Comm_rank(lugrid_tiles[g].grid.comm,&(lugrid_tiles[g].grid.iam));
            myrow = lugrid_tiles[g].grid.iam / npcol;
            mycol = lugrid_tiles[g].grid.iam % npcol;
            MPI_Comm_split(lugrid_tiles[g].grid.comm, 
                           myrow, mycol, &(lugrid_tiles[g].grid.rscp.comm));
            MPI_Comm_split(lugrid_tiles[g].grid.comm, 
                           mycol, myrow, &(lugrid_tiles[g].grid.cscp.comm));
            lugrid_tiles[g].grid.rscp.Np  = npcol;
            lugrid_tiles[g].grid.rscp.Iam = mycol;
            lugrid_tiles[g].grid.cscp.Np  = nprow;
            lugrid_tiles[g].grid.cscp.Iam = myrow;
         }
      } /* end for group g */
      if( nprocs != k )
      {
         printf("Error nprocs %d  k %d \n", nprocs, k);
         exit(-1);
      }
      free (usermap);
      solver->ML_subgroup = mygroup;
      /*
       * Fact = DOFACT Trans = NOTRANS Equil = EQUI RowPerm = LargeDiag
       * ColPerm = COLAMD ReplaceTinyPivot = REPLACE IterRefine = DOUBLE
       */
      set_default_options(optionsptr);
      optionsptr->Equil = NEQU;
      optionsptr->IterRefine = NOREFINE;
      /*
       * Future possiblities to experiment with include
       * optionsptr->RowPerm = NOROWPERM;
       * optionsptr->ColPerm = COLAMD;   (default)
       * optionsptr->ColPerm = MMD_AT_PLUS_A;
       * optionsptr->ColPerm = NATURAL;
       * ... and Equil
       */
      ScalePermstruct = ( ScalePermstruct_t *) ML_allocate( sizeof( ScalePermstruct_t));
      ScalePermstructInit(n, n, ScalePermstruct);
      LUstruct = ( LUstruct_t *) ML_allocate( sizeof( LUstruct_t) );
      LUstructInit(n, n, LUstruct);
      solver->PERMspl = (void *) ScalePermstruct;
      solver->LUspl = (void *) LUstruct;
      solver->gridtiles = lugrid_tiles;
      /* rst: mygrid is a pointer to a structure, not a structure.
       *  mygrid = ( gridinfo_t *) ML_allocate( sizeof( gridinfo_t) );
       */
   } 
   else 
   {
      /* Indicate that the factored form of A is supplied. */
      /* Reset options */
      optionsptr->Fact = FACTORED;
      optionsptr->Trans = NOTRANS;
      optionsptr->Equil = NEQU;
      optionsptr->RowPerm = MY_PERMR;
      optionsptr->ColPerm = MY_PERMC;
      optionsptr->ReplaceTinyPivot = REPLACE;
      optionsptr->IterRefine = NOREFINE;
      lugrid_tiles    = solver->gridtiles;
      mygroup = (int ) solver->ML_subgroup;
   }
   mygrid = &((lugrid_tiles[mygroup]).grid);

   /* ------------------------------------------------------------- */
   /* gather from all processors the complete right hand side       */
   /* ------------------------------------------------------------- */

   nrhs = 1;
   ML_memory_alloc((void**) &local_rhs, n*sizeof(double),"LU1" );
   for ( i = 0; i < N_local; i++ ) local_rhs[i] = rhs[i];
   i = N_local;
   ML_Comm_GappendDouble((ML_Comm *) mlcomm, local_rhs, &i, n);

   /* ------------------------------------------------------------- */
   /* perform LU decomposition and then solve                       */
   /* ------------------------------------------------------------- */

   info = flag;
   PStatInit(&stat);
   pdgssvx_ABglobal(optionsptr, A, ScalePermstruct, local_rhs, n, 
                    nrhs, mygrid, LUstruct, &berr, &stat, &info);

/*
if( iam == 0 ){
heap_info(&fragments, &total_free, &largest_free, &total_used);
printf("memory usage: fragments %d free: total %d, largest %d, total_used %d\n",
          fragments, total_free, largest_free, total_used);
}
*/

   if ( flag == 0 ) 
   {
     if ( A != NULL ) 
     {
/*
        ML_memory_free((void*)&(((NCformat *) A->Store)->rowind));
        ML_memory_free((void*)&(((NCformat *) A->Store)->colptr));
        ML_memory_free((void*)&(((NCformat *) A->Store)->nzval));
        SUPERLU_FREE( A->Store );
*/
        /* to satisfy pdgssvx_ABglobal argument check, postpone
         * ML_memory_free((void**) &A);
         */
      }
   }

   solver->reuse_flag = 1;
   PStatFree(&stat);
   if( info != 0 )
   {
      if( iam == 0 )printf("Error: ml_superlu    info = %d\n",info);
      return(-1);
   }
   /* ------------------------------------------------------------- */
   /* extract the local solution sub-vector and then clean up       */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < N_local; i++ ) x[i] = local_rhs[i+offset];
   ML_memory_free( (void **) &local_rhs );
#else
   printf("ML_SuperLU_Solve : SuperLU not used.\n");
#endif
#endif
   return 0;
}

#endif

