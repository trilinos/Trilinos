/* ************************************************************************ */
/* See the file COPYRIGHT for a complete copyright notice, contact person   */
/* and disclaimer.                                                          */
/* ************************************************************************ */

/****************************************************************************/ 
/* Declaration of the MLI_Solver functions                                  */
/****************************************************************************/ 
/* Author        : Charles Tong (LLNL)                                      */
/* Date          : December, 1999                                           */
/****************************************************************************/ 

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#ifdef ML_MPI
#include "mpi.h"
#include "ml_struct.h"
#include "ml_amg_genP.h"
#include "ml_agg_genP.h"
#include "mli_solver.h"

/****************************************************************************/ 
/* communication functions on parallel platforms                            */
/*--------------------------------------------------------------------------*/

int MLI_Irecv(void* buf, unsigned int count, int *src, int *mid,
             MPI_Comm comm, MPI_Request *request )
{
   int my_id, lsrc, retcode;

   if ( *src < 0 ) lsrc = MPI_ANY_SOURCE; else lsrc = (*src); 
   retcode = MPI_Irecv( buf, (int) count, MPI_BYTE, lsrc, *mid, comm, request);
   if ( retcode != 0 )
   {
      MPI_Comm_rank(comm, &my_id);
      printf("%d : MLI_Irecv warning : retcode = %d\n", my_id, retcode);
   }
   return 0;
}

int MLI_SIrecv(void* buf, unsigned int count, int *src, int *mid,
             MPI_Comm comm, MPI_Request *request )
{
   return 0;
}

/*--------------------------------------------------------------------------*/

int MLI_Wait(void* buf, unsigned int count, int *src, int *mid,
             MPI_Comm comm, MPI_Request *request )
{
   MPI_Status status;
   int        my_id, incount, retcode;

   retcode = MPI_Wait( request, &status );
   if ( *src < 0 ) *src = status.MPI_SOURCE; 
   if ( retcode != 0 )
   {
      MPI_Comm_rank(comm, &my_id);
      printf("%d : MLI_Wait warning : retcode = %d\n", my_id, retcode);
   }
   MPI_Get_count(&status, MPI_BYTE, &incount);
   return incount;
}

int MLI_SWait(void* buf, unsigned int count, int *src, int *mid,
             MPI_Comm comm, MPI_Request *request )
{
   return 0;
}

/*--------------------------------------------------------------------------*/

int MLI_Send(void* buf, unsigned int count, int dest, int mid, MPI_Comm comm )
{
   int my_id;
   int retcode = MPI_Send( buf, (int) count, MPI_BYTE, dest, mid, comm);
   if ( retcode != 0 )
   {
      MPI_Comm_rank(comm, &my_id);
      printf("%d : MLI_Send warning : retcode = %d\n", my_id, retcode);
   }
   return 0;
}

int MLI_SSend(void* buf, unsigned int count, int dest, int mid, MPI_Comm comm )
{
   return 0;
}

/****************************************************************************/ 
/* wrapper function for interprocessor communication for matvec and getrow  */
/*--------------------------------------------------------------------------*/

int MLI_CSRExchBdry(double *vec, void *obj)
{
   int           i, j, msgid, leng, src, dest, offset, *tempList;
   double        *dbuf;
   MLI_Context   *context;
   MLI_CSRMatrix *Amat;
   MPI_Comm    comm;
   MPI_Request *request; 

   int sendProcCnt, recvProcCnt;
   int *sendProc, *recvProc;
   int *sendLeng, *recvLeng;
   int **sendList, nRows;

   context     = (MLI_Context *) obj;
   Amat        = (MLI_CSRMatrix  *) context->Amat;
   comm        = context->comm;
   sendProcCnt = Amat->sendProcCnt;
   recvProcCnt = Amat->recvProcCnt;
   sendProc    = Amat->sendProc;
   recvProc    = Amat->recvProc;
   sendLeng    = Amat->sendLeng;
   recvLeng    = Amat->recvLeng;
   sendList    = Amat->sendList;
   nRows       = Amat->Nrows;

   request = (MPI_Request *) ML_allocate( recvProcCnt * sizeof( MPI_Request ));
   msgid = 234;
   offset = nRows;
   for ( i = 0; i < recvProcCnt; i++ )
   {
      leng = recvLeng[i] * sizeof( double );
      src  = recvProc[i];
      MLI_Irecv((void*) &(vec[offset]), leng, &src, &msgid, comm, &request[i]);
      offset += recvLeng[i];
   }
   msgid = 234;
   for ( i = 0; i < sendProcCnt; i++ )
   {
      dest = sendProc[i];
      leng = sendLeng[i] * sizeof( double );
      dbuf = (double *) ML_allocate( leng );
      tempList = sendList[i];
      for ( j = 0; j < sendLeng[i]; j++ ) {
         dbuf[j] = vec[tempList[j]];
      }
      MLI_Send((void*) dbuf, leng, dest, msgid, comm);
      if ( dbuf != NULL ) ML_free( dbuf );
   }
   offset = nRows;
   for ( i = 0; i < recvProcCnt; i++ )
   {
      leng = recvLeng[i] * sizeof( double );
      src  = recvProc[i];
      MLI_Wait((void*) &(vec[offset]), leng, &src, &msgid, comm, &request[i]);
      offset += recvLeng[i];
   }
   ML_free( request );
   return 1;
}

/****************************************************************************/ 
/* matvec function for local matrix structure MLI_CSRMatrix                 */
/*--------------------------------------------------------------------------*/

int MLI_CSRMatVec(ML_Operator *obj, int leng1, double p[], int leng2, double ap[])
{
    MLI_Context    *context;
    MLI_CSRMatrix  *Amat;

    int    i, j, length, nRows, ibeg, iend, k;
    double *dbuf, sum;
    int    *rowptr, *colnum;
    double *values;
    ML_Operator    *mat_in;

    mat_in = (ML_Operator *) obj;
    context = (MLI_Context *) ML_Get_MyMatvecData(mat_in);
    Amat    = (MLI_CSRMatrix*) context->Amat;
    nRows = Amat->Nrows;
    rowptr  = Amat->rowptr;
    colnum  = Amat->colnum;
    values  = Amat->values;

    length = nRows;
    for ( i = 0; i < Amat->recvProcCnt; i++ ) length += Amat->recvLeng[i];
    dbuf = (double *) ML_allocate( length * sizeof( double ) );
    for ( i = 0; i < nRows; i++ ) dbuf[i] = p[i];
    MLI_CSRExchBdry(dbuf, obj);
    for ( i = 0 ; i < nRows; i++ ) 
    {
       sum = 0.0;
       ibeg = rowptr[i];
       iend = rowptr[i+1];
       for ( j = ibeg; j < iend; j++ )
       { 
          k = colnum[j];
          sum += ( values[j] * dbuf[k] );
       }
       ap[i] = sum;
    }
    if ( dbuf != NULL ) ML_free( dbuf );
    return 1;
}

/****************************************************************************/
/* getrow function for local matrix structure MLI_CSRMatrix (ML compatible) */
/*--------------------------------------------------------------------------*/

int MLI_CSRGetRow(ML_Operator *obj, int N_requested_rows, int requested_rows[],
    int allocated_space, int columns[], double values[], int row_lengths[])
{
    int           i, j, ncnt, colindex, rowLeng, rowindex;
    MLI_Context   *context;
    MLI_CSRMatrix *Amat;
    int           *rowptr  = Amat->rowptr;
    int           *colInd  = Amat->colnum;
    double        *colVal  = Amat->values;
    ML_Operator   *mat_in;

    mat_in = (ML_Operator *) obj;
    context = (MLI_Context *) ML_Get_MyGetrowData(mat_in);
    Amat    = (MLI_CSRMatrix*) context->Amat;

    ncnt = 0;
    for ( i = 0; i < N_requested_rows; i++ )
    {
       rowindex = requested_rows[i];
       rowLeng = rowptr[rowindex+1] - rowptr[rowindex];
       if ( ncnt+rowLeng > allocated_space ) return 0;
       row_lengths[i] = rowLeng;
       colindex = rowptr[rowindex];
       for ( j = 0; j < rowLeng; j++ )
       {
          columns[ncnt] = colInd[colindex];
          values[ncnt++] = colVal[colindex++];
       }
    }
    return 1;
}

/****************************************************************************/
/* MLI_Solver_Create                                                        */
/*--------------------------------------------------------------------------*/

MLI_Solver *MLI_Solver_Create( MPI_Comm comm )
{
    /* create an internal ML data structure */

    MLI_Solver *solver = (MLI_Solver *) ML_allocate(sizeof(MLI_Solver));
    if ( solver == NULL ) return NULL;   

    /* fill in all other default parameters */

    solver->comm          = comm;
    solver->nlevels       = 20;   /* max number of levels */
    solver->method        = 0;    /* 0 = conjugate gradient */
    solver->pre           = -1;   /* none */
    solver->post          = -1;   /* none */
    solver->pre_sweeps    = 0;    /* 0 smoothing steps */
    solver->post_sweeps   = 2;
    solver->jacobi_wt     = 0.0;  /* default damping factor */
    solver->ml_ag         = NULL;
    solver->ml_amg        = NULL;
    solver->ag_threshold  = 0.08; /* threshold for aggregation */
    solver->ag_coarsen    = 1;    /* 1 = Uncoupled */
    solver->ag_method     = 1;    /* 1 = AMG */
    solver->contxt        = NULL; /* context for matvec */
    solver->ndiag         = 0;    /* for scaling the matrix */
    solver->diag_scale    = NULL;
    solver->nRows         = 0;
    solver->mat_ia        = NULL;
    solver->mat_ja        = NULL;
    solver->mat_a         = NULL;
    solver->rhs           = NULL;
    solver->sol           = NULL;
    solver->nPDE          = 1;
    solver->nNullVectors  = 1;
    solver->nullSpace     = NULL;
    solver->ml_ptr        = NULL;

    return solver;
}

/****************************************************************************/
/* MLI_Solver_Destroy                                                       */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Destroy( MLI_Solver *solver )
{
    int           i;
    MLI_CSRMatrix *Amat;

    if ( solver->ml_ag  != NULL ) ML_Aggregate_Destroy( &(solver->ml_ag) );
    if ( solver->ml_amg != NULL ) ML_AMG_Destroy( &(solver->ml_amg) );
    if ( solver->ml_ptr != NULL ) ML_Destroy( &(solver->ml_ptr) );
    if ( solver->contxt->partition != NULL ) ML_free( solver->contxt->partition );
    if ( solver->contxt->Amat != NULL )
    {
       Amat = (MLI_CSRMatrix *) solver->contxt->Amat;
       if ( Amat->sendProc != NULL ) ML_free(Amat->sendProc);
       if ( Amat->sendLeng != NULL ) ML_free(Amat->sendLeng);
       if ( Amat->sendList != NULL ) 
       {
          for (i = 0; i < Amat->sendProcCnt; i++ )
             if (Amat->sendList[i] != NULL) ML_free(Amat->sendList[i]);
          ML_free(Amat->sendList);
       }
       if ( Amat->recvProc != NULL ) ML_free(Amat->recvProc);
       if ( Amat->recvLeng != NULL ) ML_free(Amat->recvLeng);
       if ( Amat->map      != NULL ) ML_free(Amat->map);
       ML_free( Amat );
    }
    if ( solver->contxt != NULL ) ML_free( solver->contxt );
    if ( solver->diag_scale != NULL ) ML_free( solver->diag_scale );
    if ( solver->mat_ia    != NULL ) ML_free( solver->mat_ia );
    if ( solver->mat_ja    != NULL ) ML_free( solver->mat_ja );
    if ( solver->mat_a     != NULL ) ML_free( solver->mat_a );
    if ( solver->rhs       != NULL ) ML_free( solver->rhs   );
    if ( solver->sol       != NULL ) ML_free( solver->sol   );
    if ( solver->nullSpace != NULL ) ML_free( solver->nullSpace );
    ML_free( solver );

    return 0;
}

/****************************************************************************/
/* MLI_Solver_Setup                                                         */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Setup(MLI_Solver *solver, double *sol)
{
    int            i, my_id, nprocs, coarsest_level, level, sweeps, nlevels;
    int            *row_partition, localEqns, length, nPDE, nullDim;
    int            Nblocks, *blockList, *mat_ia, *mat_ja;
    double         *mat_a, *rhs;
    double         wght, *dble_array;
    MLI_Context    *context;
    MLI_CSRMatrix  *mli_mat;
    ML             *ml;

    /* -------------------------------------------------------- */ 
    /* create an ML object                                      */
    /* -------------------------------------------------------- */ 

    ML_Create( &(solver->ml_ptr), solver->nlevels );
    ml          = solver->ml_ptr;
    nlevels     = solver->nlevels;
    solver->sol = sol;
   
    /* -------------------------------------------------------- */ 
    /* set up the parallel environment                          */
    /* -------------------------------------------------------- */ 

    MPI_Comm_rank(solver->comm, &my_id);
    MPI_Comm_size(solver->comm, &nprocs);

    /* -------------------------------------------------------- */ 
    /* fetch the matrix row partition information and put it    */
    /* into the matrix data object (for matvec and getrow)      */
    /* -------------------------------------------------------- */ 

    localEqns = solver->nRows;
    row_partition = (int*) ML_allocate( (nprocs+1) * sizeof(int));
    for (i=0; i<nprocs; i++) row_partition[i] = 0;
    row_partition[my_id+1] = localEqns;
    MPI_Allgather(&localEqns,1,MPI_INT,&row_partition[1],1,MPI_INT,solver->comm);
    for (i=2; i<=nprocs; i++) row_partition[i] += row_partition[i-1];
    context = (MLI_Context *) ML_allocate(sizeof(MLI_Context));
    solver->contxt = context;
    context->comm = solver->comm;

    context->globalEqns = row_partition[nprocs];
    context->partition = row_partition;
    mli_mat = ( MLI_CSRMatrix * ) ML_allocate( sizeof( MLI_CSRMatrix) );
    context->Amat = mli_mat;
    mat_ia = solver->mat_ia;
    mat_ja = solver->mat_ja;
    mat_a  = solver->mat_a;
    rhs    = solver->rhs  ;
    MLI_Solver_Construct_CSRMatrix(localEqns,mat_ia,mat_ja,mat_a,mli_mat,
                            solver, context->partition,context); 

    if ( solver->ag_method == 2 ) /* scaling for smoothed aggregation */
       for ( i = 0; i < localEqns; i++ ) rhs[i] *= solver->diag_scale[i];

    /* -------------------------------------------------------- */ 
    /* set up the ML communicator information                   */
    /* -------------------------------------------------------- */ 

    ML_Set_Comm_Communicator(ml, solver->comm);
    ML_Set_Comm_MyRank(ml, my_id);
    ML_Set_Comm_Nprocs(ml, nprocs);
    ML_Set_Comm_Send(ml, MLI_Send);
    ML_Set_Comm_Recv(ml, MLI_Irecv);
    ML_Set_Comm_Wait(ml, MLI_Wait);

    /* -------------------------------------------------------- */ 
    /* set up the ML matrix information                         */
    /* -------------------------------------------------------- */ 

    ML_Init_Amatrix(ml,nlevels-1,localEqns,localEqns,(void *)context);
    MLnew_Set_Amatrix_Matvec(ml, nlevels-1, MLI_CSRMatVec);
    length = localEqns;
    for (i=0; i<mli_mat->recvProcCnt; i++) length += mli_mat->recvLeng[i];
    MLnew_Set_Amatrix_Getrow(ml,nlevels-1,MLI_CSRGetRow,MLI_CSRExchBdry,length);

    /* -------------------------------------------------------- */ 
    /* set up the mg method                                     */
    /* -------------------------------------------------------- */ 

    if ( solver->ag_method == 2 )
    {
       /* ----------------------------------------------------- */ 
       /* create an aggregate context                           */
       /* ----------------------------------------------------- */ 

       ML_Aggregate_Create(&(solver->ml_ag));
       ML_Aggregate_Set_MaxLevels( solver->ml_ag, solver->nlevels );
       ML_Aggregate_Set_Threshold( solver->ml_ag, solver->ag_threshold );
       if ( solver->ag_coarsen == 1 )
          ML_Aggregate_Set_CoarsenScheme_Uncoupled( solver->ml_ag);
       else
          ML_Aggregate_Set_CoarsenScheme_Coupled( solver->ml_ag);
       dble_array = (double*) ML_allocate( localEqns * sizeof(double));
       if ( dble_array == NULL )
       {
          printf("memory allocation problem in MLI_Solver_Setup %d.\n",localEqns);
          exit(1);
       }
       for ( i = 0; i < localEqns; i++ )
          dble_array[i] = 1.0 / solver->diag_scale[i];
       dble_array = solver->nullSpace;
       nPDE = solver->nPDE;
       nullDim = solver->nNullVectors;
       ML_Aggregate_Set_NullSpace(solver->ml_ag,nPDE,nullDim,dble_array,localEqns);
       ML_free(dble_array);

       /* ----------------------------------------------------- */ 
       /* perform aggregation                                   */
       /* ----------------------------------------------------- */ 

       coarsest_level = ML_Gen_MGHierarchy_UsingAggregation(ml, nlevels-1, 
                                        ML_DECREASING, solver->ml_ag);
    }
    else if ( solver->ag_method == 1 )
    {
       /* ----------------------------------------------------- */ 
       /* create an AMG context                                 */
       /* ----------------------------------------------------- */ 

       ML_AMG_Create(&(solver->ml_amg));
       ML_AMG_Set_MaxLevels( solver->ml_amg, solver->nlevels );
       coarsest_level = ML_Gen_MGHierarchy_UsingAMG(ml, nlevels-1,
                                        ML_DECREASING, solver->ml_amg);
    }
    if ( my_id == 0 )
       printf("ML : number of levels = %d\n", coarsest_level);

    coarsest_level = nlevels - coarsest_level;

    /* -------------------------------------------------------- */ 
    /* set up coarse grid solve                                 */
    /* -------------------------------------------------------- */ 

    if ( nlevels - 1 != coarsest_level )
    {
#ifdef SUPERLU
       ML_Gen_CoarseSolverSuperLU(ml, coarsest_level);
#else
       ML_Gen_Smoother_GaussSeidel(ml, coarsest_level, ML_PRESMOOTHER, 10, 1.0);
#endif
    } else
       ML_Gen_Smoother_OverlappedDDILUT(ml,coarsest_level,ML_PRESMOOTHER);

    /* -------------------------------------------------------- */ 
    /* set up smoother and coarse solver                        */
    /* -------------------------------------------------------- */ 

    for (level = nlevels-1; level > coarsest_level; level--) 
    {
       sweeps = solver->pre_sweeps;
       wght   = solver->jacobi_wt;
       if (wght == 0.0 ) wght = 1.0 / ml->spectral_radius[level];
       if ( my_id == 0 ) printf("level = %d, omega = %e\n", level, wght);
       switch ( solver->pre )
       {
          case 0 :
             ML_Gen_Smoother_Jacobi(ml, level, ML_PRESMOOTHER, sweeps, wght);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up presmoother : Jacobi(%d,%e)\n",sweeps,wght);
             }
             break;
          case 1 :
             ML_Gen_Smoother_GaussSeidel(ml, level, ML_PRESMOOTHER, sweeps, wght);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up presmoother : GS(%d)\n",sweeps);
             }
             break;
          case 2 :
             ML_Gen_Smoother_SymGaussSeidel(ml,level,ML_PRESMOOTHER,sweeps,wght);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up presmoother : SGS(%d,1.0)\n",sweeps);
             }
             break;
          case 3 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockJacobi(ml,level,ML_PRESMOOTHER,
                                         sweeps, wght, Nblocks, blockList);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up presmoother : VBJacobi(%d,%e)\n",sweeps,wght);
             }
             break;
          case 4 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockSymGaussSeidel(ml,level,ML_PRESMOOTHER,
                                              sweeps, wght, Nblocks, blockList);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up presmoother : VBSGS(%d)\n",sweeps);
             }
             break;
          case 5 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml,level,
                                ML_PRESMOOTHER, sweeps, wght,Nblocks, blockList);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up presmoother : VBSGSSeq(%d)\n",sweeps);
             }
             break;
          case 6 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockKrylovJacobi(ml,level,ML_PRESMOOTHER,
                                         sweeps, wght, Nblocks, blockList);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up presmoother : VBJacKrylov(%d,%e)\n",sweeps,
                        wght);
             }
             break;
          case 7 :
             ML_Gen_Smoother_OverlappedDDILUT(ml,level,ML_PRESMOOTHER);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up presmoother : Overlapped DDILUT\n");
             }
             break;
          case 8 :
             ML_Gen_Smoother_VBlockAdditiveSchwarz(ml,level,ML_PRESMOOTHER,sweeps,
                                                   0, NULL);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up presmoother : Additive Schwarz \n");
             }
             break;
          case 9 :
             ML_Gen_Smoother_VBlockMultiplicativeSchwarz(ml,level,ML_PRESMOOTHER,
                                                         sweeps, 0, NULL);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up presmoother : Multiplicative Schwarz \n");
             }
             break;
          default :
             printf("Setting up presmoother : invalid smoother \n");
             exit(1);
             break;
       }

       sweeps = solver->post_sweeps;
       switch ( solver->post )
       {
          case 0 :
             ML_Gen_Smoother_Jacobi(ml, level, ML_POSTSMOOTHER, sweeps, wght);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up postmoother : Jacobi(%d,%e)\n",sweeps,wght);
             }
             break;
          case 1 :
             ML_Gen_Smoother_GaussSeidel(ml, level, ML_POSTSMOOTHER, sweeps,wght);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up postsmoother : GS(%d)\n",sweeps);
             }
             break;
          case 2 :
             ML_Gen_Smoother_SymGaussSeidel(ml,level,ML_POSTSMOOTHER,sweeps,wght);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up postsmoother : SGS(%d,1.0)\n",sweeps);
             }
             break;
          case 3 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockJacobi(ml,level,ML_POSTSMOOTHER,
                                         sweeps, wght, Nblocks, blockList);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up postsmoother : VBJacobi(%d,%e)\n",sweeps,wght);
             }
             break;
          case 4 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockSymGaussSeidel(ml,level,ML_POSTSMOOTHER,
                                          sweeps, wght, Nblocks, blockList);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up postsmoother : VBSGS(%d)\n",sweeps);
             }
             break;
          case 5 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml,level,
                            ML_PRESMOOTHER, sweeps, wght, Nblocks, blockList);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up postsmoother : VBSGSSeq(%d)\n",sweeps);
             }
             break;
          case 6 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockKrylovJacobi(ml,level,ML_POSTSMOOTHER,
                                         sweeps, wght, Nblocks, blockList);
             if ( my_id == 0 && level == nlevels-1)
             {
                printf("Setting up postsmoother : VBJacKrylov(%d,%e)\n",sweeps,
                        wght);
             }
             break;
          default :
             printf("Setting up postsmoother : invalid smoother \n");
             exit(1);
             break;
       }
    }

    ML_Gen_Solver(ml, ML_MGV, nlevels-1, coarsest_level);
   
    return 0;
}

/****************************************************************************/
/* MLI_Solver_Solve                                                         */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Solve( MLI_Solver *solver)
{
    double      *diag_scale, *sol, *rhs, elapsed_time;
    ML          *ml = solver->ml_ptr;
    int         i, leng, level = ml->ML_num_levels - 1;
    ML_Operator *Amat = &(ml->Amat[level]);
    ML_Krylov   *ml_kry;

    ml_kry = ML_Krylov_Create(ml->comm);
    ML_Krylov_Set_Method(ml_kry, solver->method);
    ML_Krylov_Set_Amatrix(ml_kry, Amat);
    ML_Krylov_Set_Precon(ml_kry, ml);
    ML_Krylov_Set_PreconFunc(ml_kry, ML_AMGVSolve_Wrapper);

    leng = Amat->outvec_leng;
    diag_scale = solver->diag_scale;

    /********* another choice : diagonal preconditioner 
    ML_Krylov_Set_Precon(ml_kry, ml_kry);
    ML_Krylov_Set_PreconFunc(ml_kry, ML_DiagScale_Wrapper);
    diag = (double *) ML_allocate(leng * sizeof(double));
    for ( i = 0; i < leng; i++ )
       diag[i] = diag_scale[i] * diag_scale[i];
    ML_Krylov_Set_Diagonal(ml_kry, leng, diag);
    ML_free( diag );
    *******************************************************/

    rhs = solver->rhs;
    sol = solver->sol;
    elapsed_time = GetClock();
    ML_Krylov_Solve(ml_kry, leng, rhs, sol);
    elapsed_time = GetClock() - elapsed_time;
    if ( ml->comm->ML_mypid == 0 )
       printf("====> Solve Time = %e\n", elapsed_time);
    ML_Krylov_Destroy(&ml_kry);
  
    if ( solver->ag_method == 2 )
       for ( i = 0; i < leng; i++ ) sol[i] *= diag_scale[i];
    return 0;
}

/****************************************************************************/
/* MLI_Solver_Set_MLNumLevels                                               */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Set_MLNumLevels(MLI_Solver *solver,int nlevels)
{
    if ( nlevels < 0 )
    {
       printf("MLI_Solver_Set_MLNumLevels error : nlevels set to 1.\n");
       solver->nlevels = 1;
    } 
    else
    {
       solver->nlevels = nlevels;
    } 
    return( 0 );
}

/****************************************************************************/
/* MLI_Solver_Set_KrylovMethod                                           */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Set_KrylovMethod(MLI_Solver *solver,int method)
{
    if ( method < 0 || method > 2 )
    {
       printf("MLI_Solver_Set_KrylovMethod error : reset to CG.\n");
       solver->method = 0;
    } 
    else
    {
       solver->method = method;
    } 
    return( 0 );
}

/****************************************************************************/
/* MLI_Solver_Set_StrongThreshold                                           */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Set_StrongThreshold(MLI_Solver *solver,double strong_threshold)
{
    if ( strong_threshold < 0.0 )
    {
       printf("MLI_Solver_Set_StrongThreshold error : reset to 0.\n");
       solver->ag_threshold = 0.0;
    } 
    else
    {
       solver->ag_threshold = strong_threshold;
    } 
    return( 0 );
}

/****************************************************************************/
/* MLI_Solver_Set_NumPreSmoothings                                          */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Set_NumPreSmoothings( MLI_Solver *solver, int num_sweeps  )
{
    if ( num_sweeps < 0 )
    {
       printf("MLI_Solver_Set_NumPreSmoothings error : reset to 0.\n");
       solver->pre_sweeps = 0;
    } 
    else
    {
       solver->pre_sweeps = num_sweeps;
    } 
    return( 0 );
}

/****************************************************************************/
/* MLI_Solver_Set_PostSmoothings                                          */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Set_NumPostSmoothings( MLI_Solver *solver, int num_sweeps  )
{
    if ( num_sweeps < 0 )
    {
       printf("MLI_Solver_Set_NumPostSmoothings error : reset to 0.\n");
       solver->post_sweeps = 0;
    } 
    else
    {
       solver->post_sweeps = num_sweeps;
    } 
    return( 0 );
}

/****************************************************************************/
/* MLI_Solver_Set_PreSmoother                                               */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Set_PreSmoother( MLI_Solver *solver, int smoother_type  )
{
    if ( smoother_type < 0 || smoother_type > 7 )
    {
       printf("MLI_Solver_Set_PreSmoother error : set to Jacobi.\n");
       solver->pre = 0;
    } 
    else
    {
       solver->pre = smoother_type;
    } 
    return( 0 );
}

/****************************************************************************/
/* MLI_Solver_Set_PostSmoother                                              */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Set_PostSmoother( MLI_Solver *solver, int smoother_type  )
{
    if ( smoother_type < 0 || smoother_type > 7 )
    {
       printf("MLI_Solver_Set_PostSmoother error : set to Jacobi.\n");
       solver->post = 0;
    } 
    else
    {
       solver->post = smoother_type;
    } 
    return( 0 );
}

/****************************************************************************/
/* MLI_Solver_Set_DampingFactor                                             */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Set_DampingFactor( MLI_Solver *solver, double factor  )
{
    if ( factor < 0.0 || factor > 2.0 )
    {
       printf("MLI_Solver_Set_DampingFactor error : set to 0.5.\n");
       solver->jacobi_wt = 0.5;
    } 
    else
    {
       solver->jacobi_wt = factor;
    } 
    return( 0 );
}

/****************************************************************************/
/* MLI_Solver_Set_CoarsenScheme                                             */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Set_CoarsenScheme( MLI_Solver *solver, int scheme  )
{
    if ( scheme != 1 && scheme != 2 )
    {
       printf("MLI_Solver_Set_CoarsenScheme error : scheme set to Uncoupled.\n");
       solver->ag_coarsen = 1;
    } 
    else
    {
       solver->ag_coarsen = scheme;
    } 
    return( 0 );
}

/****************************************************************************/
/* MLI_Solver_Set_MGMethod                                                  */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Set_MGMethod( MLI_Solver *solver, int method  )
{
    if ( method != 1 && method != 2 )
    {
       printf("MLI_Solver_Set_MGMethod error : method set to AMG.\n");
       solver->ag_method  = 1;
    } 
    else
    {
       solver->ag_method = method;
    } 
    return( 0 );
}

/****************************************************************************/
/* MLI_Solver_Construct_CSRMatrix                                           */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Construct_CSRMatrix(int nrows, int *mat_ia, int *mat_ja,
                                   double *mat_a, MLI_CSRMatrix *mli_mat,
                                   MLI_Solver *solver, int *partition,
                                   MLI_Context *obj) 
{
    int         i, j, index, my_id, nprocs, msgid, *tempCnt;
    int         sendProcCnt, *sendLeng, *sendProc, **sendList;
    int         recvProcCnt, *recvLeng, *recvProc, ncnt, nnz;
    int         rowLeng, *colInd, startRow, endRow, localEqns;
    int         *diagSize, *offdiagSize, externLeng, *externList;
    int         num_bdry;
    double      *colVal, *diag_scale;
    MPI_Status  status;
    MPI_Comm    comm;

    /* -------------------------------------------------------- */
    /* get machine information and local matrix information     */
    /* -------------------------------------------------------- */
    
    comm = solver->comm;
    MPI_Comm_rank(comm, &my_id);
    MPI_Comm_size(comm, &nprocs);

    startRow  = partition[my_id];
    endRow    = partition[my_id+1] - 1;
    localEqns = endRow - startRow + 1;

    /* -------------------------------------------------------- */
    /* probe A to find out about diagonal and off-diagonal      */
    /* block information                                        */
    /* -------------------------------------------------------- */

    diagSize    = (int*) ML_allocate( sizeof(int) * localEqns );
    offdiagSize = (int*) ML_allocate( sizeof(int) * localEqns );
    num_bdry = 0;
    for ( i = startRow; i <= endRow; i++ )
    {
       diagSize[i-startRow] = offdiagSize[i-startRow] = 0;
       colInd = &(mat_ja[mat_ia[i-startRow]]);
       colVal = &(mat_a[mat_ia[i-startRow]]);
       rowLeng = mat_ia[i-startRow+1] - mat_ia[i-startRow];
       for (j = 0; j < rowLeng; j++)
       {
          if ( colInd[j] < startRow || colInd[j] > endRow )
          {
             offdiagSize[i-startRow]++;
          }
          else
          {
             diagSize[i-startRow]++;
          }
       }
       if ( diagSize[i-startRow] + offdiagSize[i-startRow] == 1 ) 
          num_bdry++;
    }

    /* -------------------------------------------------------- */
    /* construct external node list in global eqn numbers       */
    /* -------------------------------------------------------- */

    externLeng = 0;
    for ( i = 0; i < localEqns; i++ ) externLeng += offdiagSize[i];
    if ( externLeng > 0 )
         externList = (int *) ML_allocate( sizeof(int) * externLeng);
    else externList = NULL;
    externLeng = 0;
    for ( i = startRow; i <= endRow; i++ )
    {
       colInd = &(mat_ja[mat_ia[i-startRow]]);
       colVal = &(mat_a[mat_ia[i-startRow]]);
       rowLeng = mat_ia[i-startRow+1] - mat_ia[i-startRow];
       for (j = 0; j < rowLeng; j++)
       {
          if ( colInd[j] < startRow || colInd[j] > endRow )
             externList[externLeng++] = colInd[j];
       }
    }
    if ( externLeng > 0 ) ML_sort( externLeng, externList );
    ncnt = 1;
    for ( i = 1; i < externLeng; i++ )
    {
       if ( externList[i] != externList[ncnt-1] ) 
          externList[ncnt++] = externList[i];
    }
    if ( externLeng > 0 ) externLeng = ncnt;
    else                  externLeng = 0;

    /* -------------------------------------------------------- */
    /* allocate the CSR matrix                                  */
    /* -------------------------------------------------------- */ 

    diag_scale = (double *) ML_allocate((nrows+externLeng)*sizeof(double)); 
    nnz = 0; 
    for (i = 0; i < localEqns; i++) nnz += diagSize[i] + offdiagSize[i]; 
    ML_free( diagSize );
    ML_free( offdiagSize );

    /* -------------------------------------------------------- */ 
    /* put the matrix data in the CSR matrix                    */
    /* -------------------------------------------------------- */ 

    for ( i = startRow; i <= endRow; i++ )
    {
       colInd = &(mat_ja[mat_ia[i-startRow]]);
       colVal = &(mat_a[mat_ia[i-startRow]]);
       rowLeng = mat_ia[i-startRow+1] - mat_ia[i-startRow];
       for (j = 0; j < rowLeng; j++)
       {
          index = colInd[j];
          if ( index < startRow || index > endRow )
          {
             colInd[j] = ML_sorted_search(index,externLeng,externList );
             colInd[j] += localEqns;
          }
          else
          {
             if ( colInd[j] == i ) {
                diag_scale[i-startRow] = colVal[j];
             }
             colInd[j] = index - startRow;
          }
       }
    }
   
    /* -------------------------------------------------------- */ 
    /* initialize the MLI_CSRMatrix data structure              */
    /* -------------------------------------------------------- */ 

    mli_mat->Nrows       = localEqns;
    mli_mat->rowptr      = mat_ia;
    mli_mat->colnum      = mat_ja;
    mli_mat->values      = mat_a;
    mli_mat->sendProcCnt = 0;
    mli_mat->recvProcCnt = 0;
    mli_mat->sendLeng    = NULL;
    mli_mat->recvLeng    = NULL;
    mli_mat->sendProc    = NULL;
    mli_mat->recvProc    = NULL;
    mli_mat->sendList    = NULL;
    mli_mat->map         = externList;
 
    /* -------------------------------------------------------- */ 
    /* form the remote portion of the matrix                    */
    /* -------------------------------------------------------- */ 

    if ( nprocs > 1 ) 
    {
       /* ----------------------------------------------------- */ 
       /* count number of elements to be received from each     */
       /* remote processor (assume sequential mapping)          */
       /* ----------------------------------------------------- */ 

       tempCnt = (int *) ML_allocate( sizeof(int) * nprocs );
       for ( i = 0; i < nprocs; i++ ) tempCnt[i] = 0;
       for ( i = 0; i < externLeng; i++ )
       {
          for ( j = 0; j < nprocs; j++ )
          {
             if ( externList[i] >= partition[j] && 
                  externList[i] < partition[j+1] )
             {
                tempCnt[j]++;
                break;
             }
          }
       }

       /* ----------------------------------------------------- */ 
       /* compile a list processors data is to be received from */
       /* ----------------------------------------------------- */ 

       recvProcCnt = 0;
       for ( i = 0; i < nprocs; i++ )
          if ( tempCnt[i] > 0 ) recvProcCnt++;
       recvLeng = (int*) ML_allocate( sizeof(int) * recvProcCnt );
       recvProc = (int*) ML_allocate( sizeof(int) * recvProcCnt );
       recvProcCnt = 0;
       for ( i = 0; i < nprocs; i++ )
       {
          if ( tempCnt[i] > 0 ) 
          {
             recvProc[recvProcCnt]   = i;
             recvLeng[recvProcCnt++] = tempCnt[i];
          }
       }

       /* ----------------------------------------------------- */ 
       /* each processor has to find out how many processors it */
       /* has to send data to                                   */
       /* ----------------------------------------------------- */ 

       sendLeng = (int *) ML_allocate( nprocs * sizeof(int) );
       for ( i = 0; i < nprocs; i++ ) tempCnt[i] = 0;
       for ( i = 0; i < recvProcCnt; i++ ) tempCnt[recvProc[i]] = 1;
       MPI_Allreduce(tempCnt, sendLeng, nprocs, MPI_INT, MPI_SUM, comm );
       sendProcCnt = sendLeng[my_id];
       ML_free( sendLeng );
       if ( sendProcCnt > 0 )
       {
          sendLeng = (int *)  ML_allocate( sendProcCnt * sizeof(int) );
          sendProc = (int *)  ML_allocate( sendProcCnt * sizeof(int) );
          sendList = (int **) ML_allocate( sendProcCnt * sizeof(int*) );
       }
       else 
       {
          sendLeng = sendProc = NULL;
          sendList = NULL;
       }

       /* ----------------------------------------------------- */ 
       /* each processor sends to all processors it expects to  */
       /* receive data about the lengths of data expected       */
       /* ----------------------------------------------------- */ 

       msgid = 539;
       for ( i = 0; i < recvProcCnt; i++ ) 
       {
          MPI_Send((void*) &recvLeng[i],1,MPI_INT,recvProc[i],msgid,comm);
       }
       for ( i = 0; i < sendProcCnt; i++ ) 
       {
          MPI_Recv((void*) &sendLeng[i],1,MPI_INT,MPI_ANY_SOURCE,msgid,
                   comm,&status);
          sendProc[i] = status.MPI_SOURCE;
          sendList[i] = (int *) ML_allocate( sendLeng[i] * sizeof(int) );
          if ( sendList[i] == NULL ) 
             printf("allocate problem %d \n", sendLeng[i]);
       }

       /* ----------------------------------------------------- */ 
       /* each processor sends to all processors it expects to  */
       /* receive data about the equation numbers               */
       /* ----------------------------------------------------- */ 

       for ( i = 0; i < nprocs; i++ ) tempCnt[i] = 0; 
       ncnt = 1;
       for ( i = 0; i < externLeng; i++ ) 
       {
          if ( externList[i] >= partition[ncnt] )
          {
             tempCnt[ncnt-1] = i;
             i--;
             ncnt++;
          }
       }    
       for ( i = ncnt-1; i < nprocs; i++ ) tempCnt[i] = externLeng; 

       /* ----------------------------------------------------- */ 
       /* send the global equation numbers                      */
       /* ----------------------------------------------------- */ 

       msgid = 540;
       for ( i = 0; i < recvProcCnt; i++ ) 
       {
          if ( recvProc[i] == 0 ) j = 0;
          else                    j = tempCnt[recvProc[i]-1];
          rowLeng = recvLeng[i];
          MPI_Send((void*) &externList[j],rowLeng,MPI_INT,recvProc[i],
                    msgid,comm);
       }
       for ( i = 0; i < sendProcCnt; i++ ) 
       {
          rowLeng = sendLeng[i];
          MPI_Recv((void*)sendList[i],rowLeng,MPI_INT,sendProc[i],
                   msgid,comm,&status);
       }

       /* ----------------------------------------------------- */ 
       /* convert the send list from global to local numbers    */
       /* ----------------------------------------------------- */ 

       for ( i = 0; i < sendProcCnt; i++ )
       { 
          for ( j = 0; j < sendLeng[i]; j++ )
          {
             index = sendList[i][j] - startRow;
             if ( index < 0 || index >= localEqns )
             {
                printf("%d: Construct MLI matrix Error - index",my_id);
                printf(" out of range%d\n", index);
             }
             sendList[i][j] = index;
          }
       }

       /* ----------------------------------------------------- */ 
       /* convert the send list from global to local numbers    */
       /* ----------------------------------------------------- */ 

       mli_mat->sendProcCnt = sendProcCnt;
       mli_mat->recvProcCnt = recvProcCnt;
       mli_mat->sendLeng    = sendLeng;
       mli_mat->recvLeng    = recvLeng;
       mli_mat->sendProc    = sendProc;
       mli_mat->recvProc    = recvProc;
       mli_mat->sendList    = sendList;

       /* ----------------------------------------------------- */ 
       /* clean up                                              */
       /* ----------------------------------------------------- */ 

       ML_free( tempCnt );
    }

    /* -------------------------------------------------------- */ 
    /* finally scale the matrix                                 */
    /* -------------------------------------------------------- */ 

    MLI_CSRExchBdry(diag_scale, obj);
    if ( solver->ag_method == 2 )
    {
       for ( i = 0; i < nrows+externLeng; i++ )
       {
          if ( diag_scale[i] != 0.0 ) 
             diag_scale[i] = 1.0 / sqrt(diag_scale[i]);
       }
       for ( i = 0; i < nrows; i++ )
       {
          for ( j = mat_ia[i]; j < mat_ia[i+1]; j++ )
             mat_a[j] = mat_a[j] * diag_scale[i] * diag_scale[mat_ja[j]];
       }
       solver->diag_scale = diag_scale;
    }
    else
    {
       ML_free( diag_scale );
       solver->diag_scale = NULL;
    }
    return 0;
}

/****************************************************************************/
/* reading a matrix from a file in ija format (first row : nrows, nnz)      */
/* (read by a single processor)                                             */
/*--------------------------------------------------------------------------*/

void MLI_Solver_Read_IJAFromFile(double **val, int **ia, int **ja, int *N,
                                double **rhs, char *matfile, char *rhsfile)
{
   int    i, Nrows, nnz, icount, rowindex, colindex, curr_row;
   int    k, m, *mat_ia, *mat_ja, ncnt, rnum;
   double dtemp, *mat_a, value, *rhs_local;
   FILE   *fp;

   /* ----------------------------------------------------------- */
   /* read matrix file                                            */
   /* ----------------------------------------------------------- */

   printf("Reading matrix file = %s \n", matfile );
   fp = fopen( matfile, "r" );
   if ( fp == NULL ) {
      printf("Error : file open error (filename=%s).\n", matfile);
      exit(1);
   }
   fscanf(fp, "%d %d", &Nrows, &nnz);
   if ( Nrows <= 0 || nnz <= 0 ) {
      printf("Error : nrows,nnz = %d %d\n", Nrows, nnz);
      exit(1);
   }
   mat_ia = (int *)    ML_allocate((Nrows+1) * sizeof(int));
   mat_ja = (int *)    ML_allocate(nnz * sizeof(int));
   mat_a  = (double *) ML_allocate(nnz * sizeof(double));
   mat_ia[0] = 0;

   curr_row = 0;
   icount   = 0;
   for ( i = 0; i < nnz; i++ ) {
      fscanf(fp, "%d %d %lg", &rowindex, &colindex, &value);
      rowindex--;
      colindex--;
      if ( rowindex != curr_row ) mat_ia[++curr_row] = icount;
      if ( rowindex < 0 || rowindex >= Nrows )
         printf("Error reading row %d (curr_row = %d)\n", rowindex, curr_row);
      if ( colindex < 0 || colindex >= Nrows )
         printf("Error reading col %d (rowindex = %d)\n", colindex, rowindex);
      if ( rowindex == colindex && value == 0.0 ) 
         printf("Warning : zero diagonal at row %d\n", rowindex);
      /*
      if ( value != 0.0 ) {
      */      
         mat_ja[icount] = colindex;
         mat_a[icount++]  = value;
      /*
      }
      */      
      if ( i % 100000 == 0 ) printf("processing data %d\n", i);
   }
   fclose(fp);
   for ( i = curr_row+1; i <= Nrows; i++ ) mat_ia[i] = icount;
   (*val) = mat_a;
   (*ia)  = mat_ia;
   (*ja)  = mat_ja;
   (*N) = Nrows;
   printf("matrix has %6d rows and %7d nonzeros\n", Nrows, mat_ia[Nrows]);

   /* ----------------------------------------------------------- */
   /* read rhs file                                               */
   /* ----------------------------------------------------------- */

   printf("reading rhs file = %s \n", rhsfile );
   fp = fopen( rhsfile, "r" );
   if ( fp == NULL ) {
      printf("No rhs file found.\n");
      rhs_local = (double *) ML_allocate(Nrows * sizeof(double));
      for ( k = 0; k < Nrows; k++ ) {
         rhs_local[k] = drand48();
      }
      (*rhs) = rhs_local;
   }
   else
   {
      fscanf(fp, "%d", &ncnt);
      if ( ncnt <= 0 || ncnt != Nrows) {
         printf("Error : nrows = %d \n", ncnt);
         exit(1);
      }
      fflush(stdout);
      rhs_local = (double *) ML_allocate(Nrows * sizeof(double));
      m = 0;
      for ( k = 0; k < ncnt; k++ ) {
         fscanf(fp, "%d %lg", &rnum, &dtemp);
         rhs_local[rnum-1] = dtemp; m++;
      }
      fflush(stdout);
      ncnt = m;
      fclose(fp);
      (*rhs) = rhs_local;
      printf("reading rhs done \n");
   }

   printf("returning \n");
}

/****************************************************************************/
/* reading a matrix from a file in ija format (first row : nrows, nnz)      */
/* (read by a single processor)                                             */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Get_IJAFromFile(MLI_Solver *solver, char *matfile, char *rhsfile)
{
   int      i, j, chunksize, mybegin, myend, my_id, nprocs, nnz;
   int      blksize=3, local_N, ncnt;
   int      mat_N, *mat_ia, *mat_ja, *ia, *ja;
   double   *mat_a, *mat_rhs, *val, *rhs;
   MPI_Comm comm;

   comm = solver->comm;

   MPI_Comm_rank(comm, &my_id);
   MPI_Comm_size(comm, &nprocs);

   if ( my_id == 0 ) {
      MLI_Solver_Read_IJAFromFile(&mat_a, &mat_ia, &mat_ja, &mat_N,
                                  &mat_rhs, "matrix.data", "rhs.data");
      nnz = mat_ia[mat_N];
      MPI_Bcast(&mat_N, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&nnz,   1, MPI_INT, 0, MPI_COMM_WORLD);

      MPI_Bcast(mat_ia,  mat_N+1, MPI_INT,    0, MPI_COMM_WORLD);
      MPI_Bcast(mat_ja,  nnz,     MPI_INT,    0, MPI_COMM_WORLD);
      MPI_Bcast(mat_a,   nnz,     MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(mat_rhs, mat_N,   MPI_DOUBLE, 0, MPI_COMM_WORLD);

   } else {
      MPI_Bcast(&mat_N, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&nnz,   1, MPI_INT, 0, MPI_COMM_WORLD);
      mat_ia  = (int    *) ML_allocate( (mat_N + 1) * sizeof( int ) );
      mat_ja  = (int    *) ML_allocate( nnz * sizeof( int ) );
      mat_a   = (double *) ML_allocate( nnz * sizeof( double ) );
      mat_rhs = (double *) ML_allocate( mat_N * sizeof( double ) );

      MPI_Bcast(mat_ia,  mat_N+1, MPI_INT,    0, MPI_COMM_WORLD);
      MPI_Bcast(mat_ja,  nnz,     MPI_INT,    0, MPI_COMM_WORLD);
      MPI_Bcast(mat_a,   nnz,     MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(mat_rhs, mat_N,   MPI_DOUBLE, 0, MPI_COMM_WORLD);
   }

   chunksize = mat_N / blksize;
   if ( chunksize * blksize != mat_N )
   {
      printf("Cannot put into matrix blocks with block size 3\n");
      exit(1);
   }
   chunksize = chunksize / nprocs;
   mybegin = chunksize * my_id * blksize;
   myend   = chunksize * (my_id + 1) * blksize - 1;
   if ( my_id == nprocs-1 ) myend = mat_N - 1;
   printf("Processor %d : begin/end = %d %d\n", my_id, mybegin, myend);
   fflush(stdout);

   local_N = myend - mybegin + 1;
   nnz = 0;
   for ( i = mybegin; i <= myend; i++ ) nnz += (mat_ia[i+1] - mat_ia[i]);
   ia  = (int *) ML_allocate( (local_N+1) * sizeof(int));
   ja  = (int *) ML_allocate( nnz * sizeof(int));
   val = (double *) ML_allocate( nnz * sizeof(double));
   rhs = (double *) ML_allocate( local_N * sizeof(double));
   ia[0] = 0;
   ncnt = 0;
   for ( i = mybegin; i <= myend; i++ ) 
   {
      for ( j = mat_ia[i]; j < mat_ia[i+1]; j++ ) 
      {
         ja[ncnt] = mat_ja[j];
         val[ncnt++] = mat_a[j];
      }
      ia[i-mybegin+1] = ncnt;
      rhs[i-mybegin] = mat_rhs[i];
   }
   solver->nRows = local_N;
   solver->mat_ia = ia;
   solver->mat_ja = ja;
   solver->mat_a  = val;
   solver->rhs    = rhs;
   ML_free( mat_ia);
   ML_free( mat_ja);
   ML_free( mat_a);
   return local_N;
} 

/****************************************************************************/
/* reading a number of vectors corresponding to the rigid body mode         */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Get_NullSpaceFromFile(MLI_Solver *solver, char *rbmfile)
{
   int      i, k, chunksize, mybegin, myend, my_id, nprocs;
   int      blksize=3, local_N, ncnt, global_N, nullDimension;
   double   dtemp, *rbm_global, *rbm;
   MPI_Comm comm;
   FILE     *fp;

   comm = solver->comm;

   MPI_Comm_rank(comm, &my_id);
   MPI_Comm_size(comm, &nprocs);

   local_N = solver->nRows;
   MPI_Allreduce(&local_N, &global_N, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   if ( my_id == 0 ) {
      fp = fopen( rbmfile, "r" );
      if ( fp == NULL ) {
         printf("No rbm file found.\n");
         if ( global_N == 0 )
         {
            if ( my_id == 0 ) 
               printf("Error : should call MatrixRead before RBMRead\n");
            exit(1);
         }
         rbm_global = (double *) ML_allocate(global_N * sizeof(double));
         for ( k = 0; k < global_N; k++ ) rbm_global[k] = 1.0;
         nullDimension = 1;
      }
      else
      {
         fscanf(fp, "%d", &nullDimension);
         fscanf(fp, "%d", &ncnt);
         if ( ncnt != global_N ) {
            printf("Error in reading RBM : dimension not matched.\n");
            exit(1);
         }
         rbm_global = (double *) ML_allocate(global_N*nullDimension*sizeof(double));
         for ( k = 0; k < ncnt*nullDimension; k++ ) {
            fscanf(fp, "%lg", &dtemp);
            rbm_global[k] = dtemp;
         }
      }
   }
 
   chunksize = global_N / blksize;
   chunksize = chunksize / nprocs;
   mybegin = chunksize * my_id * blksize;
   myend   = chunksize * (my_id + 1) * blksize - 1;
   if ( my_id == nprocs-1 ) myend = global_N - 1;
   if ( myend - mybegin + 1 != local_N )
   {
      printf("Error : RBM and matrix local dimensions do not match.\n");
      printf("        info = %d %d\n", myend-mybegin+1,local_N);
      exit(1);
   }

   MPI_Bcast(&nullDimension, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if ( my_id != 0 )
      rbm_global = (double *) ML_allocate(global_N*nullDimension*sizeof(double));

   MPI_Bcast(rbm_global,global_N*nullDimension,MPI_DOUBLE,0,MPI_COMM_WORLD);

   rbm = (double *) ML_allocate(local_N * nullDimension * sizeof(double));
   ncnt = 0;
   for ( k = 0; k < nullDimension; k++ ) 
      for ( i = mybegin; i <= myend; i++ ) 
         rbm[ncnt++] = rbm_global[k*global_N+i];

   if ( nullDimension == 1 ) solver->nPDE = 1;
   else                      solver->nPDE = 3;
   solver->nNullVectors = nullDimension;
   solver->nullSpace = rbm;
   ML_free( rbm_global);
   return mybegin;
} 

/****************************************************************************/
/* MLI_Solver_SetupDD                                                       */
/*--------------------------------------------------------------------------*/

int MLI_Solver_SetupDD(MLI_Solver *solver,int startRow, int Nrows,
                     int *mat_ia,int *mat_ja, double *mat_a, double *b, 
                     double *x)
{
    int            i, my_id, nprocs, coarsest_level, level, sweeps, nlevels;
    int            *row_partition, localEqns, length;
    int            Nblocks, *blockList, one=1;
    double         wght;
    MLI_Context    *context;
    MLI_CSRMatrix  *mli_mat;

    /* -------------------------------------------------------- */ 
    /* fetch the ML pointer                                     */
    /* -------------------------------------------------------- */ 

    ML      *ml   = solver->ml_ptr;
    nlevels       = solver->nlevels;
   
    /* -------------------------------------------------------- */ 
    /* set up the serial environment                            */
    /* -------------------------------------------------------- */ 

    MPI_Comm_rank(solver->comm, &my_id);
    MPI_Comm_size(solver->comm, &nprocs);

    /* -------------------------------------------------------- */ 
    /* fetch the matrix row partition information and put it    */
    /* into the matrix data object (for matvec and getrow)      */
    /* -------------------------------------------------------- */ 

    localEqns = Nrows;
    row_partition = (int*) ML_allocate( (nprocs+1) * sizeof(int));
    for (i=0; i<nprocs; i++) row_partition[i] = 0;
    row_partition[my_id+1] = Nrows;
    MPI_Allgather(&Nrows,1,MPI_INT,&row_partition[1],1,MPI_INT,solver->comm);
    for (i=2; i<=nprocs; i++) row_partition[i] += row_partition[i-1];
    context = (MLI_Context *) ML_allocate(sizeof(MLI_Context));
    solver->contxt = context;
    context->comm = (USR_COMM) one;

    context->globalEqns = Nrows;
    context->partition = NULL;
    mli_mat = ( MLI_CSRMatrix * ) ML_allocate( sizeof( MLI_CSRMatrix) );
    context->Amat = mli_mat;
    MLI_Solver_Construct_LocalCSRMatrix(Nrows,mat_ia,mat_ja,mat_a,
                         mli_mat, solver, context->partition,context); 
    for ( i = 0; i < Nrows; i++ ) b[i] *= solver->diag_scale[i];

    /* -------------------------------------------------------- */ 
    /* set up the ML communicator information                   */
    /* -------------------------------------------------------- */ 

    ML_Set_Comm_Communicator(ml, (USR_COMM) one);
    ML_Set_Comm_MyRank(ml, 0);
    ML_Set_Comm_Nprocs(ml, 1);
    ML_Set_Comm_Send(ml, MLI_SSend);
    ML_Set_Comm_Recv(ml, MLI_SIrecv);
    ML_Set_Comm_Wait(ml, MLI_SWait);

    /* -------------------------------------------------------- */ 
    /* set up the ML matrix information                         */
    /* -------------------------------------------------------- */ 

    ML_Init_Amatrix(ml,nlevels-1,localEqns,localEqns,(void *)context);
    MLnew_Set_Amatrix_Matvec(ml, nlevels-1, MLI_CSRMatVec);
    length = localEqns;
    for (i=0; i<mli_mat->recvProcCnt; i++) length += mli_mat->recvLeng[i];
    MLnew_Set_Amatrix_Getrow(ml,nlevels-1,MLI_CSRGetRow,MLI_CSRExchBdry,length);

    /* -------------------------------------------------------- */ 
    /* create an aggregate context                              */
    /* -------------------------------------------------------- */ 

    ML_Aggregate_Create(&(solver->ml_ag));
    solver->ml_ag->max_levels = solver->nlevels;
    ML_Aggregate_Set_Threshold( solver->ml_ag, solver->ag_threshold );

    /* -------------------------------------------------------- */ 
    /* perform aggregation                                      */
    /* -------------------------------------------------------- */ 

    coarsest_level = ML_Gen_MGHierarchy_UsingAggregation(ml, nlevels-1, 
                                        ML_DECREASING, solver->ml_ag);
    if ( my_id == 0 )
       printf("ML : number of levels = %d\n", coarsest_level);

    coarsest_level = nlevels - coarsest_level;

    /* -------------------------------------------------------- */ 
    /* set up smoother and coarse solver                        */
    /* -------------------------------------------------------- */ 

    for (level = nlevels-1; level > coarsest_level; level--) 
    {
       sweeps = solver->pre_sweeps;
       wght   = solver->jacobi_wt;
       switch ( solver->pre )
       {
          case 0 :
             ML_Gen_Smoother_Jacobi(ml, level, ML_PRESMOOTHER, sweeps, wght);
             break;
          case 1 :
             ML_Gen_Smoother_GaussSeidel(ml,level,ML_PRESMOOTHER,sweeps,wght);
             break;
          case 2 :
             ML_Gen_Smoother_SymGaussSeidel(ml,level,ML_PRESMOOTHER,sweeps,wght);
             break;
          case 3 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockJacobi(ml,level,ML_PRESMOOTHER,
                                         sweeps, wght, Nblocks, blockList);
             break;
          case 4 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockSymGaussSeidel(ml,level,ML_PRESMOOTHER,
                                              sweeps, wght, Nblocks, blockList);
             break;
          case 5 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml,level,
                                ML_PRESMOOTHER, sweeps, wght,Nblocks, blockList);
             break;
          case 6 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockKrylovJacobi(ml,level,ML_PRESMOOTHER,
                                         sweeps, wght, Nblocks, blockList);
             break;
          case 7 :
             ML_Gen_Smoother_OverlappedDDILUT(ml,level,ML_PRESMOOTHER);
             break;
          case 8 :
             ML_Gen_Smoother_VBlockAdditiveSchwarz(ml,level,ML_PRESMOOTHER,sweeps,
                                                   0, NULL);
             break;
          case 9 :
             ML_Gen_Smoother_VBlockMultiplicativeSchwarz(ml,level,ML_PRESMOOTHER,
                                                         sweeps, 0, NULL);
             break;
          default :
             printf("Setting up presmoother : invalid smoother \n");
             exit(1);
             break;
       }

       sweeps = solver->post_sweeps;
       switch ( solver->post )
       {
          case 0 :
             ML_Gen_Smoother_Jacobi(ml, level, ML_POSTSMOOTHER, sweeps, wght);
             break;
          case 1 :
             ML_Gen_Smoother_GaussSeidel(ml,level,ML_POSTSMOOTHER,sweeps,wght);
             break;
          case 2 :
             ML_Gen_Smoother_SymGaussSeidel(ml,level,ML_POSTSMOOTHER,sweeps,wght);
             break;
          case 3 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockJacobi(ml,level,ML_POSTSMOOTHER,
                                         sweeps, wght, Nblocks, blockList);
             break;
          case 4 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockSymGaussSeidel(ml,level,ML_POSTSMOOTHER,
                                          sweeps, wght, Nblocks, blockList);
             break;
          case 5 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockSymGaussSeidelSequential(ml,level,
                            ML_POSTSMOOTHER, sweeps, wght,Nblocks, blockList);
             break;
          case 6 :
             Nblocks = ML_Aggregate_Get_AggrCount( solver->ml_ag, level );
             ML_Aggregate_Get_AggrMap( solver->ml_ag, level, &blockList );
             ML_Gen_Smoother_VBlockKrylovJacobi(ml,level,ML_POSTSMOOTHER,
                                         sweeps, wght, Nblocks, blockList);
             break;
          case 7 :
             ML_Gen_Smoother_OverlappedDDILUT(ml,level,ML_POSTSMOOTHER);
             break;
          default :
             printf("Setting up postsmoother : invalid smoother \n");
             exit(1);
             break;
       }
    }

    ML_Gen_CoarseSolverSuperLU(ml, coarsest_level);
    /*ML_Gen_Smoother_GaussSeidel(ml, coarsest_level, ML_PRESMOOTHER, 100);*/

    ML_Gen_Solver(ml, ML_MGV, nlevels-1, coarsest_level);
   
    return 0;
}

/****************************************************************************/
/* MLI_Solver_Construct_CSRMatrix                                           */
/*--------------------------------------------------------------------------*/

int MLI_Solver_Construct_LocalCSRMatrix(int nrows, int *mat_ia, int *mat_ja,
                                   double *mat_a, MLI_CSRMatrix *mli_mat,
                                   MLI_Solver *solver, int *partition,
                                   MLI_Context *obj) 
{
    int         i, j, index, ncnt, nnz;
    int         rowLeng, *colInd, startRow, endRow, localEqns;
    int         *diagSize, num_bdry;
    double      *colVal, *diag_scale;

    /* -------------------------------------------------------- */
    /* get matrix information                                   */
    /* -------------------------------------------------------- */
    
    startRow  = 0;
    endRow    = nrows - 1;
    localEqns = endRow - startRow + 1;

    /* -------------------------------------------------------- */
    /* probe A to find out about diagonal and off-diagonal      */
    /* block information                                        */
    /* -------------------------------------------------------- */

    diagSize    = (int*) ML_allocate( sizeof(int) * localEqns );
    num_bdry = 0;
    for ( i = startRow; i <= endRow; i++ )
    {
       diagSize[i-startRow] = 0;
       colInd = &(mat_ja[mat_ia[i-startRow]]);
       colVal = &(mat_a[mat_ia[i-startRow]]);
       rowLeng = mat_ia[i-startRow+1] - mat_ia[i-startRow];
       for (j = 0; j < rowLeng; j++)
       {
          if ( colInd[j] >= startRow && colInd[j] <= endRow )
          {
             diagSize[i-startRow]++;
          }
       }
       if ( diagSize[i-startRow] == 1 ) num_bdry++;
    }

    /* -------------------------------------------------------- */
    /* allocate the CSR matrix                                  */
    /* -------------------------------------------------------- */ 

    diag_scale = (double *) ML_allocate(nrows*sizeof(double)); 
    nnz = 0; 
    for (i = 0; i < localEqns; i++) nnz += diagSize[i]; 
    ML_free( diagSize );

    /* -------------------------------------------------------- */ 
    /* put the matrix data in the CSR matrix                    */
    /* -------------------------------------------------------- */ 

    ncnt = 0;
    mat_ia[0] = ncnt;
    for ( i = startRow; i <= endRow; i++ )
    {
       colInd = &(mat_ja[mat_ia[i-startRow]]);
       colVal = &(mat_a[mat_ia[i-startRow]]);
       rowLeng = mat_ia[i-startRow+1] - mat_ia[i-startRow];
       for (j = 0; j < rowLeng; j++)
       {
          index = colInd[j];
          if ( index >= startRow && index <= endRow )
          {
             if ( index == i ) {
                diag_scale[i-startRow] = colVal[j];
             }
             mat_ja[ncnt]  = index - startRow;
             mat_a[ncnt++] = colVal[j];
          }
       }
       mat_ia[i-startRow+1] = ncnt;
    }
   
   
    /* -------------------------------------------------------- */ 
    /* initialize the MLI_CSRMatrix data structure              */
    /* -------------------------------------------------------- */ 

    mli_mat->Nrows       = localEqns;
    mli_mat->rowptr      = mat_ia;
    mli_mat->colnum      = mat_ja;
    mli_mat->values      = mat_a;
    mli_mat->sendProcCnt = 0;
    mli_mat->recvProcCnt = 0;
    mli_mat->sendLeng    = NULL;
    mli_mat->recvLeng    = NULL;
    mli_mat->sendProc    = NULL;
    mli_mat->recvProc    = NULL;
    mli_mat->sendList    = NULL;
    mli_mat->map         = NULL;
 
    /* -------------------------------------------------------- */ 
    /* finally scale the matrix                                 */
    /* -------------------------------------------------------- */ 

    for ( i = 0; i < nrows; i++ )
    {
       if ( diag_scale[i] != 0.0 ) 
           diag_scale[i] = 1.0 / sqrt(diag_scale[i]);
    }
    for ( i = 0; i < nrows; i++ )
    {
       for ( j = mat_ia[i]; j < mat_ia[i+1]; j++ )
          mat_a[j] = mat_a[j] * diag_scale[i] * diag_scale[mat_ja[j]];
    }
    solver->diag_scale = diag_scale;
    
    return 0;
}
#endif

