/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/*       User Interface Functions                                            */
/* ************************************************************************* */
/* ************************************************************************* */

#include <math.h>
#include "ml_struct.h"
#include "ml_aztec_utils.h"
#include "ml_agg_genP.h"
#include "ml_lapack.h"
#include "ml_smoother.h"
#ifdef ML_MPI
#include "mpi.h"
#endif

#ifdef SUPERLU
extern int SuperLU_Solve(void *, int, double *, int, double *);
extern int SuperLU_SolveLocal(void *, double *, double *);
#elif DSUPERLU
extern int SuperLU_Solve(void *, int, double *, int, double *);
extern int SuperLU_SolveLocal(void *, double *, double *);
#endif


/* ************************************************************************* */
/* create and initialize a ML object                                         */
/* ------------------------------------------------------------------------- */

int ML_Create(ML **ml_ptr, int Nlevels)
{
   int             i, length;
   double          *max_eigen;
   ML_Operator     *Amat, *Rmat, *Pmat;
   ML_Smoother     *pre_smoother, *post_smoother;
   ML_CSolve       *csolve;
   ML_Grid         *Grid;
   ML_BdryPts      *BCs;
   ML_Mapper       *eqn2grid, *grid2eqn;
   ML_DVector      *Amat_Normalization;
   ML_1Level       *SingleLevel;
   char            str[80];

#ifdef ML_TIMING
   struct ML_Timing *timing;
#endif
 
   ML_memory_alloc( (void**) ml_ptr, sizeof(ML), "MLM" );

   (*ml_ptr)->output_level    = 10;
   (*ml_ptr)->res_output_freq = 1;
   (*ml_ptr)->tolerance       = 1.e-8;
   (*ml_ptr)->max_iterations  = 1000;

   ML_Comm_Create( &((*ml_ptr)->comm) );
   global_comm = (*ml_ptr)->comm;

   ML_memory_alloc((void**) &pre_smoother, sizeof(ML_Smoother)*Nlevels,"MS1");
   ML_memory_alloc((void**) &post_smoother,sizeof(ML_Smoother)*Nlevels,"MS2");
   ML_memory_alloc((void**) &csolve       ,sizeof(ML_CSolve  )*Nlevels,"MCS");
   ML_memory_alloc((void**) &Grid         ,sizeof(ML_Grid    )*Nlevels,"MGD");
   ML_memory_alloc((void**) &BCs         ,sizeof(ML_BdryPts  )*Nlevels,"MBC");
   ML_memory_alloc((void**) &eqn2grid    ,sizeof(ML_Mapper   )*Nlevels,"MM1");
   ML_memory_alloc((void**) &grid2eqn    ,sizeof(ML_Mapper   )*Nlevels,"MM2");
   ML_memory_alloc((void**) &SingleLevel ,sizeof(ML_1Level   )*Nlevels,"MSL");
   ML_memory_alloc((void**) &Amat         ,sizeof(ML_Operator)*Nlevels,"MAM");
   ML_memory_alloc((void**) &Rmat         ,sizeof(ML_Operator)*Nlevels,"MRM");
   ML_memory_alloc((void**) &Pmat         ,sizeof(ML_Operator)*Nlevels,"MPM");
   ML_memory_alloc((void**) &max_eigen    ,sizeof(ML_Operator)*Nlevels,"MQM");
   length = sizeof(ML_DVector) * Nlevels;
   for ( i = 0; i < Nlevels; i++ ) max_eigen[i] = 0.0;
   ML_memory_alloc((void**)&Amat_Normalization, length, "MAN");

   (*ml_ptr)->ML_num_levels      = Nlevels;
   (*ml_ptr)->pre_smoother       = pre_smoother;
   (*ml_ptr)->post_smoother      = post_smoother;
   (*ml_ptr)->csolve             = csolve;
   (*ml_ptr)->Amat               = Amat;
   (*ml_ptr)->Grid               = Grid;
   (*ml_ptr)->BCs                = BCs;
   (*ml_ptr)->eqn2grid           = eqn2grid;
   (*ml_ptr)->grid2eqn           = grid2eqn;
   (*ml_ptr)->SingleLevel        = SingleLevel;
   (*ml_ptr)->Rmat               = Rmat;
   (*ml_ptr)->Pmat               = Pmat;
   (*ml_ptr)->spectral_radius = max_eigen;
   (*ml_ptr)->Amat_Normalization = Amat_Normalization ;
   (*ml_ptr)->timing             = NULL;

#ifdef ML_TIMING
   ML_memory_alloc((void**) &timing, sizeof(struct ML_Timing),"MT");
   timing->precond_apply_time = 0.;
   timing->total_build_time   = 0.;
   (*ml_ptr)->timing = timing;
#endif 

   for (i = 0; i < Nlevels; i++) 
   {
      ML_Operator_Init(&(Amat[i]), (*ml_ptr)->comm);
      ML_Operator_Set_1Levels(&(Amat[i]), &SingleLevel[i], &SingleLevel[i]);
      ML_Operator_Set_BdryPts(&(Amat[i]), &BCs[i]);
      ML_Operator_Init(&(Rmat[i]), (*ml_ptr)->comm);
      ML_Operator_Set_1Levels(&(Rmat[i]), &SingleLevel[i], NULL);
      ML_Operator_Set_BdryPts(&(Rmat[i]), &BCs[i]);
      ML_Operator_Init(&(Pmat[i]), (*ml_ptr)->comm);
      ML_Operator_Set_1Levels(&(Pmat[i]), &SingleLevel[i], NULL);
      ML_Operator_Set_BdryPts(&(Pmat[i]), NULL);

      (SingleLevel[i]).comm = (ML_Comm *) (*ml_ptr)->comm;
      SingleLevel[i].Amat          = &Amat[i];
      SingleLevel[i].Rmat          = &Rmat[i];
      SingleLevel[i].Pmat          = &Pmat[i];
      SingleLevel[i].BCs           = &BCs[i];
      SingleLevel[i].eqn2grid      = &eqn2grid[i];
      SingleLevel[i].grid2eqn      = &grid2eqn[i];
      SingleLevel[i].Grid          = &Grid[i];
      SingleLevel[i].pre_smoother  = &pre_smoother[i];
      SingleLevel[i].post_smoother = &post_smoother[i];
      SingleLevel[i].csolve        = &csolve[i];
      SingleLevel[i].Amat_Normalization = &Amat_Normalization[i];
      ML_DVector_Init( &Amat_Normalization[i] );
      SingleLevel[i].levelnum      = i;

      ML_Mapper_Init( &(eqn2grid[i]) );
      ML_Mapper_Init( &(grid2eqn[i]) );
      ML_Grid_Init( &(Grid[i]) );
      ML_BdryPts_Init( &(BCs[i]) );

      ML_Smoother_Init( &(pre_smoother[i]), &(SingleLevel[i]) );
      ML_Smoother_Init( &(post_smoother[i]), &(SingleLevel[i]) );

      ML_CSolve_Init( &(csolve[i]) );
      ML_CSolve_Set_1Level( &(csolve[i]), &(SingleLevel[i]) );
      sprintf(str,"Amat_%d",i); ML_Operator_Set_Label( &(Amat[i]),str);
      sprintf(str,"Rmat_%d",i); ML_Operator_Set_Label( &(Rmat[i]),str);
      sprintf(str,"Pmat_%d",i); ML_Operator_Set_Label( &(Pmat[i]),str);
      sprintf(str,"PreS_%d",i); ML_Smoother_Set_Label( &(pre_smoother[i]),str);
      sprintf(str,"PostS_%d",i);ML_Smoother_Set_Label( &(post_smoother[i]),str);
      sprintf(str,"Solve_%d",i);ML_CSolve_Set_Label(&(csolve[i]),str);
  }
  ML_random_init();
  return 0;
}

/* ************************************************************************* */
/* destroy an ML object                                                      */
/* ------------------------------------------------------------------------- */

int ML_Destroy(ML **ml_ptr)
{
   int i;
   ML  *ml;
 
   ml = (*ml_ptr);

#ifdef ML_TIMING
   if ( ml->output_level != 0 ) ML_Print_Timing(ml);
#endif

   for (i = 0; i < ml->ML_num_levels; i++) 
   {
      ML_Operator_Clean(&(ml->Amat[i]));
      ML_Operator_Clean(&(ml->Rmat[i]));
      ML_Operator_Clean(&(ml->Pmat[i]));
      ML_Grid_Clean(&(ml->Grid[i]));
      ML_BdryPts_Clean(&(ml->BCs[i]));
      ML_DVector_Clean( &(ml->Amat_Normalization[i]) );
      ML_Smoother_Clean(&(ml->pre_smoother[i]));
      ML_Smoother_Clean(&(ml->post_smoother[i]));
      ML_CSolve_Clean(&(ml->csolve[i]));
   }

   ML_memory_free( (void**) &(ml->csolve[0].func ) );
   ML_memory_free( (void**) &(ml->pre_smoother) );
   ML_memory_free( (void**) &(ml->post_smoother) );
   ML_memory_free( (void**) &(ml->csolve) );
   ML_memory_free( (void**) &(ml->Amat) );
   ML_memory_free( (void**) &(ml->Rmat) );
   ML_memory_free( (void**) &(ml->Pmat) );
   ML_memory_free( (void**) &(ml->Amat_Normalization) );
   ML_memory_free( (void**) &(ml->Grid) );
   ML_memory_free( (void**) &(ml->BCs) );
   ML_memory_free( (void**) &(ml->eqn2grid) );
   ML_memory_free( (void**) &(ml->grid2eqn) );
   ML_memory_free( (void**) &(ml->SingleLevel) );
   ML_memory_free( (void**) &(ml->spectral_radius) );
   if (ml->timing != NULL) ML_memory_free( (void**) &(ml->timing) );
   ML_Comm_Destroy( &(ml->comm) );
   ML_memory_free( (void**) &(ml) );
   (*ml_ptr) = NULL;
#ifdef ML_DEBUG
   ML_memory_inquire();
#endif
   return 0;
}

/* ************************************************************************* */
/* set debug level                                                           */
/* ------------------------------------------------------------------------- */

int ML_Set_OutputLevel(ML *ml, int output_level)
{
  ml->output_level  = output_level;
  return(1);
}

/* ------------------------------------------------------------------------- */

int ML_Set_ResidualOutputFrequency(ML *ml, int output_freq)
{
  ml->res_output_freq = output_freq;
  return(1);
}

/* ------------------------------------------------------------------------- */

int ML_Set_Tolerance(ML *ml, double tolerance)
{
  ml->tolerance = tolerance;
  return(1);
}

/* ------------------------------------------------------------------------- */

int ML_Set_MaxIterations(ML *ml, int iterations)
{
  ml->max_iterations = iterations;
  return(1);
}

/* ************************************************************************* */
/* functions to initialize the communicator                                  */
/* ------------------------------------------------------------------------- */

int ML_Init_Comm(ML *ml)
{
   return(ML_Comm_Create( &(ml->comm) ));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Comm_MyRank(ML *ml, int myrank)
{
   return (ML_Comm_Set_Mypid( ml->comm, myrank));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Comm_Nprocs(ML *ml, int nprocs)
{
   return (ML_Comm_Set_Nprocs( ml->comm, nprocs));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Comm_Communicator(ML *ml, USR_COMM com)
{
   return (ML_Comm_Set_UsrComm(ml->comm, com));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Comm_Send(ML *ml, int (*send)())
{
   return ( ML_Comm_Set_SendFcn (ml->comm, send ) );
}

/* ------------------------------------------------------------------------- */

int ML_Set_Comm_Recv(ML *ml, int (*recv)())
{
   return ( ML_Comm_Set_RecvFcn (ml->comm, recv ) );
}

/* ------------------------------------------------------------------------- */

int ML_Set_Comm_Wait(ML *ml, int (*wait)())
{
   return ( ML_Comm_Set_WaitFcn (ml->comm, wait ) );
}

/* ------------------------------------------------------------------------- */

int ML_Set_Comm(ML *ml, ML_Comm *comm )
{
   int i;

   ml->comm = comm;
   for (i = 0; i < ml->ML_num_levels; i++) 
      ml->SingleLevel[i].comm = comm;
   return 0;
}

/* ************************************************************************* */
/* functions to initialize grids                                             */
/* ------------------------------------------------------------------------- */

int ML_Init_Grid( ML *ml, int level, void *grid)
{
   ML_Grid_Set_Grid( ml->SingleLevel[level].Grid, grid);
   ML_Grid_Create_GridFunc( ml->SingleLevel[level].Grid);
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GridFunc( ML *ml, int level, ML_GridFunc *gf)
{
   if ( ml->SingleLevel[level].Grid->gridfcn != NULL &&
        ml->SingleLevel[level].Grid->gf_SetOrLoad == 2) 
      ML_GridFunc_Destroy(&(ml->SingleLevel[level].Grid->gridfcn));
      
   ml->SingleLevel[level].Grid->gridfcn = gf;
   ml->SingleLevel[level].Grid->gf_SetOrLoad = 1; 
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid( ML *ml, int level, void *grid, ML_GridFunc *gf)
{
   ml->SingleLevel[level].Grid->Grid = grid;
   if ( ml->SingleLevel[level].Grid->gridfcn != NULL &&
        ml->SingleLevel[level].Grid->gf_SetOrLoad == 2) 
      ML_GridFunc_Destroy(&(ml->SingleLevel[level].Grid->gridfcn));
   ml->SingleLevel[level].Grid->gridfcn = gf;
   ml->SingleLevel[level].Grid->gf_SetOrLoad = 1; 
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_MaxVertPerElmnt(ML *ml, int nl, int nvert)
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_MaxVertPerElmnt( gf, nvert );
   return 0;
}
   
/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetDimension(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetDimension(gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetNVert(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetNVert( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetNElmnt(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetNElmnts( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntNVert(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntNVert( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntVertList(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntVertList( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntGlobalNum(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntGlobalNum( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetVertGlobalNum(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetVertGlobalNum( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetVertCoordinate(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetVertCoordinate( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_ComputeBasisCoef(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_ComputeBasisCoef( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntVolume(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntVolumes( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntMatrix(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntMatrix( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntNullSpace(ML *ml, int nl, int (*func)())
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntNullSpace( gf, func );
   return 0;
}

/* ************************************************************************* */
/* ************************************************************************* */
/* functions to initialize the Amatrix                                       */
/* ------------------------------------------------------------------------- */

int ML_Init_Amatrix(ML *ml, int level, int ilen, int olen, void *data)
{
   ML_Operator_Set_1Levels(&(ml->Amat[level]),&(ml->SingleLevel[level]),
			  &(ml->SingleLevel[level]));
   ML_Operator_Set_ApplyFuncData(&(ml->Amat[level]), ilen, olen, ML_EMPTY,
                             data, olen, NULL, 0);
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Amatrix_Matvec(ML *ml, int level, 
                      int (*matvec)(void *, int, double *, int, double *)) 
{
   ML_Operator *matrix;
   matrix = &(ml->Amat[level]);

   return(ML_Operator_Set_ApplyFunc(matrix,ML_EXTERNAL,matvec));
}
int ML_Get_Amatrix(ML *ml, int level, ML_Operator **matrix)
{
   *matrix = &(ml->Amat[level]);
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Amatrix_Diag(ML *ml, int nl, int size, double diagonal[])
{
   return(ML_Operator_Set_Diag(&(ml->Amat[nl]), size, diagonal) );
}

/* ------------------------------------------------------------------------- */

int ML_Set_Amatrix_Getrow(ML *ml, int nl, 
        int (*getrow)(void *, int , int* , int , int*, double* , int*),
	int (*comm  )(double *vec, void *data), int comm_vec_leng )
{
   ML_Operator *Amat;
   int         Nghost;

   Amat = &(ml->Amat[nl]);

   if (comm != NULL) {
      Nghost = comm_vec_leng - Amat->invec_leng;
      if (Nghost < 0) {
         printf("ML_Set_Amatrix_Getrow: comm_vec_leng is less than the\n");
         printf("                       matrix's invec_length\n");
         exit(1);
      }
      ML_CommInfoOP_Generate( &(Amat->getrow->pre_comm), comm, Amat->data, 
			      ml->comm, Amat->invec_leng, Nghost);
   }
   else {
      if ((ml->comm->ML_nprocs > 1) & (ml->comm->ML_mypid == 0)) {
         printf("Warning: No communication information given to ");
         printf("ML_Set_Amatrix_Getrow\n");
      }
      ML_CommInfoOP_Set_neighbors(&(Amat->getrow->pre_comm), 0,
                               NULL, ML_OVERWRITE, NULL, 0);

   }

   return(ML_Operator_Set_Getrow(Amat, ML_EXTERNAL, Amat->outvec_leng, getrow));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Amatrix_GetrowNeighbors(ML *ml, int level, int N_neigh, 
                                   int *neigh_list)
{
   ML_CommInfoOP_Set_neighbors(&(ml->SingleLevel[level].Amat->getrow->pre_comm),
			       N_neigh, neigh_list, ML_OVERWRITE, NULL, 0);
   return(1);
}

/* ------------------------------------------------------------------------- */

int ML_Set_Amatrix_GetrowCommInfo(ML *ml, int level, int neighbor, 
           int N_rcv, int *rcv_list, int N_send, int *send_list)
{
   ML_CommInfoOP_Set_exch_info(ml->SingleLevel[level].Amat->getrow->pre_comm, 
			       neighbor, N_rcv, rcv_list, N_send, send_list);
   return(1);
}

/* ------------------------------------------------------------------------- */

int ML_Set_Amatrix_NormalizationFactors(ML *ml,int level,int leng,
                                        double *data)
{
   ML_DVector_LoadData(ml->SingleLevel[level].Amat_Normalization,leng,data);
   return 0;
}

/* ************************************************************************* */
/* ************************************************************************* */
/* functions to initialize the restrictor                                    */
/* ------------------------------------------------------------------------- */

int ML_Init_Restrictor(ML *ml, int level, int level2, int ilen, int olen, 
                       void *data)
{
   ML_Operator_Set_1Levels(&(ml->Rmat[level]),&(ml->SingleLevel[level]),
                          &(ml->SingleLevel[level2]));
   ML_Operator_Set_ApplyFuncData(&(ml->Rmat[level]), ilen, olen, ML_EMPTY,
                             data, olen, NULL, 0);
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Restrictor_Matvec( ML *ml , int from_level, 
	int (*func) (void *, int, double *, int, double *))
{
   ML_Operator *matrix;
   matrix = &(ml->Rmat[from_level]);

   return(ML_Operator_Set_ApplyFunc(matrix,ML_EXTERNAL,func));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Restrictor_Getrow(ML *ml, int nl, 
        int (*getrow)(void *, int , int* , int , int*, double* , int*),
        int (*comm  )(double *vec, void *data), int comm_vec_leng )
{
   ML_Operator *Rmat;
   int         Nghost;

   Rmat = &(ml->Rmat[nl]);

   if (comm != NULL) {
      Nghost = comm_vec_leng - Rmat->invec_leng;
      if (Nghost < 0) {
         printf("ML_Set_Restrictor_Getrow: comm_vec_leng is less than the\n");
         printf("                       matrix's invec_length\n");
         exit(1);
      }
      ML_CommInfoOP_Generate( &(Rmat->getrow->pre_comm), comm, Rmat->data,
                              ml->comm, Rmat->invec_leng, Nghost);
   }
   else {
      if ((ml->comm->ML_nprocs > 1) & (ml->comm->ML_mypid == 0)) {
         printf("Warning: No communication information given to ");
         printf("ML_Set_Restrictor_Getrow\n");
      }
   }

   return(ML_Operator_Set_Getrow(Rmat, ML_EXTERNAL, Rmat->outvec_leng, getrow));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Restrictor_GetrowNeighbors(ML *ml, int level, int N_neigh, 
                                      int *neigh_list)
{
   ML_CommInfoOP_Set_neighbors(&(ml->SingleLevel[level].Rmat->getrow->pre_comm),
			       N_neigh, neigh_list, ML_OVERWRITE, NULL, 0);
   return(1);
}

/* ------------------------------------------------------------------------- */

int ML_Set_Restrictor_GetrowCommInfo(ML *ml, int level, int neighbor, 
            int N_rcv, int *rcv_list, int N_send, int *send_list)
{
   ML_CommInfoOP_Set_exch_info(ml->SingleLevel[level].Rmat->getrow->pre_comm, 
			       neighbor, N_rcv, rcv_list, N_send, send_list);
   return(1);
}

/* ************************************************************************* */
/* ************************************************************************* */
/* functions to initialize the prolongator                                   */
/* ------------------------------------------------------------------------- */

int ML_Init_Prolongator(ML *ml, int level, int level2, int ilen, int olen, 
                        void *data)
{
   ML_Operator_Set_1Levels(&(ml->Pmat[level]),&(ml->SingleLevel[level]),
                          &(ml->SingleLevel[level2]));
   ML_Operator_Set_ApplyFuncData(&(ml->Pmat[level]), ilen, olen, ML_EMPTY,
                             data, olen, NULL, 0);
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Prolongator_Matvec( ML *ml , int to_level, 
	int (*func) (void *, int, double *, int, double *))
{
   ML_Operator *matrix;
   matrix = &(ml->Pmat[to_level]);

   return(ML_Operator_Set_ApplyFunc(matrix,ML_EXTERNAL, func));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Prolongator_Getrow(ML *ml, int nl, 
        int (*getrow)(void *, int , int* , int , int*, double* , int*),
        int (*comm  )(double *vec, void *data), int comm_vec_leng )
{
   ML_Operator *Pmat;
   int         Nghost;

   Pmat = &(ml->Pmat[nl]);

   if (comm != NULL) {
      Nghost = comm_vec_leng - Pmat->invec_leng;
      if (Nghost < 0) {
         printf("ML_Set_Prolongator_Getrow: comm_vec_leng is less than \n");
         printf("                           the matrix's invec_length\n");
         exit(1);
      }
      ML_CommInfoOP_Generate( &(Pmat->getrow->pre_comm), comm, Pmat->data,
                              ml->comm, Pmat->invec_leng, Nghost);
   }
   else {
      if ((ml->comm->ML_nprocs > 1) & (ml->comm->ML_mypid == 0)) {
         printf("Warning: No communication information given to ");
         printf("ML_Set_Prolongator_Getrow\n");
      }
   }

   return(ML_Operator_Set_Getrow(Pmat, ML_EXTERNAL, Pmat->outvec_leng, getrow));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Prolongator_GetrowNeighbors(ML *ml, int level, int N_neigh, 
                                       int *neigh_list)
{
   ML_CommInfoOP_Set_neighbors(&(ml->SingleLevel[level].Pmat->getrow->pre_comm),
				N_neigh, neigh_list, ML_OVERWRITE, NULL, 0);
   return(1);
}

/* ------------------------------------------------------------------------- */

int ML_Set_Prolongator_GetrowCommInfo(ML *ml, int level, int neighbor, 
                int N_rcv, int *rcv_list, int N_send, int *send_list)
{
   ML_CommInfoOP_Set_exch_info(ml->SingleLevel[level].Pmat->getrow->pre_comm, 
			       neighbor, N_rcv, rcv_list, N_send, send_list);
   return(1);
}

/* ************************************************************************* */
/* ************************************************************************* */
/* functions to initialize the smoother                                      */
/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------- */
/* set the user-defined smoother                                             */
/* ------------------------------------------------------------------------- */

int ML_Set_Smoother( ML *ml , int nl , int pre_or_post, void *data, 
                     int (*func)(void *, int, double *, int, double *) ) 
{
   if (nl == ML_ALL_LEVELS) { 
      printf("ML_Set_Smoother: ML_ALL_LEVELS not allowed\n");
      return 1;
   }
   if (nl < 0) {
      printf("ML_Set_Smoother: cannot set smoother on level %d\n",nl);
      return 1;
   }
   if (pre_or_post == ML_PRESMOOTH) {
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_EXTERNAL, data,
                              NULL, func, 1,(double) ML_DEFAULT));
   }
   else if (pre_or_post == ML_POSTSMOOTH)
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_EXTERNAL, data,
                              NULL, func, 1, (double) ML_DEFAULT));
   else if (pre_or_post == ML_BOTH)  {
      ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_EXTERNAL, data,
                              NULL, func, 1,(double) ML_DEFAULT);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_EXTERNAL, data,
                              NULL, func, 1, (double) ML_DEFAULT));
   }
   else return(pr_error("ML_Set_Smoother: unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the damped Jacobi smoother                                       */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_Jacobi( ML *ml , int nl, int pre_or_post, int ntimes,
        double omega) 
{
   int (*fun)(void *, int, double *, int, double *);
   int start_level, end_level, i, status;

   if (nl == ML_ALL_LEVELS) { start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_Smoother_Jacobi: cannot set smoother on level %d\n",start_level);
      return 1;
   }

   fun = ML_Smoother_Jacobi;
   if (omega == ML_DEFAULT) omega = .5;

   if (pre_or_post == ML_PRESMOOTH) {
      for (i = start_level; i <= end_level; i++)
       status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL, NULL,
                             fun, NULL, ntimes, omega);
   }
   else if (pre_or_post == ML_POSTSMOOTH)
      for (i = start_level; i <= end_level; i++)
       status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL, NULL,
                             fun, NULL, ntimes,omega);
   else if (pre_or_post == ML_BOTH) {
      for (i = start_level; i <= end_level; i++) {
        status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL, NULL,
                             fun, NULL, ntimes, omega);
        status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL, NULL,
                             fun, NULL, ntimes,omega);
      }
   }
   else return(pr_error("ML_Gen_Smoother_Jacobi: unknown pre_or_post choice\n"));

   return(status);
}

/* ------------------------------------------------------------------------- */
/* generate the Gauss Seidel smoother (SOR)                                  */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_GaussSeidel( ML *ml , int nl, int pre_or_post, int ntimes,
                                double omega)
{
   int (*fun)(void *, int, double *, int, double *);
   int start_level, end_level, i, status;

   if (nl == ML_ALL_LEVELS) { start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_Smoother_GaussSeidel: cannot set smoother on level %d\n",start_level);
      return 1;
   }

   fun = ML_Smoother_GaussSeidel;

   if (pre_or_post == ML_PRESMOOTH) {
      for (i = start_level; i <= end_level; i++)
             status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL, NULL,
                                      fun, NULL, ntimes, omega);
   }
   else if (pre_or_post == ML_POSTSMOOTH) {
      for (i = start_level; i <= end_level; i++)
             status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL, NULL,
                             fun, NULL, ntimes, omega);
   }
   else if (pre_or_post == ML_BOTH) {
      for (i = start_level; i <= end_level; i++) {
             status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL, NULL,
                                      fun, NULL, ntimes, omega);
             status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL, NULL,
                             fun, NULL, ntimes, omega);
      }
   }
   else return(pr_error("ML_Gen_Gauss-Seidel: unknown pre_or_post choice\n"));
   return(status);
}

/* ------------------------------------------------------------------------- */
/* generate the symmetric Gauss Seidel smoother                              */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_SymGaussSeidel( ML *ml , int nl, int pre_or_post, 
                                   int ntimes, double omega)
{
   int         start_level, end_level, i, j, status, Nrows, count;
   int         *bindx;
   double      *nums, **sgs_nums = NULL, *num2;
   double      *val = NULL, temp_omega;
   ML_Operator *Amat;
   int         (*fun)(void *, int, double *, int, double *);
   void        (*fun2)(void *) = NULL;
   struct ML_CSR_MSRdata *ptr = NULL;

   if (nl == ML_ALL_LEVELS) { start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_Smoother_SymGaussSeidel: cannot set smoother on level %d\n",start_level);
      return 1;
   }
	
   for (i = start_level; i <= end_level; i++) {

      if (omega == ML_DEFAULT) omega = 1.;
      fun  = ML_Smoother_SGS;
      Amat = &(ml->Amat[i]);

      if (Amat->getrow->external == MSR_getrows){
         ptr   = (struct ML_CSR_MSRdata *) Amat->data;
         val   = ptr->values;
         bindx = ptr->columns;
      }

#ifdef AZTEC
      else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif
      if (val != NULL) {
         fun = ML_Smoother_MSR_SGS;
         sgs_nums = (double **) ML_allocate( sizeof(double)*2);
         Nrows    = Amat->getrow->Nrows;
         nums     = (double *) ML_allocate( Nrows * sizeof(double));
         num2     = (double *) ML_allocate( Nrows * sizeof(double));
         sgs_nums[0] = nums;
         sgs_nums[1] = num2;
         for (j = 0; j < Nrows; j++) { 
            count = 0;
         /*
         for (j = bindx[j]; j < bindx[j+1]; j++) if (bindx[j] >= Nrows) count++;
         */

            if (bindx[j] != bindx[j+1]) 
               temp_omega = omega*(1.0 - .5*((double) count) / 
                               ((double) (bindx[j+1]-bindx[j])));
            else temp_omega = 1.;

            num2[j] = 1. - temp_omega;
            if (val[j] != 0.0) nums[j] = temp_omega/val[j];
            else { nums[j] = 0.0; num2[j] = 1.; }
         }
         fun2 = ML_Smoother_Clean_MSR_GS;
      }

      if (pre_or_post == ML_PRESMOOTH) {
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL, 
                                  sgs_nums, fun, NULL, ntimes, omega);
	 ml->pre_smoother[i].data_destroy = fun2;
      }
      else if (pre_or_post == ML_POSTSMOOTH) {
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL, 
                                  sgs_nums, fun, NULL, ntimes, omega);
	 ml->post_smoother[i].data_destroy = fun2;
      }
      else if (pre_or_post == ML_BOTH) {
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL, 
                                  sgs_nums, fun, NULL, ntimes, omega);
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL, 
                                  sgs_nums, fun, NULL, ntimes, omega);
	 ml->post_smoother[i].data_destroy = fun2;
      }
      else return(pr_error("Print unknown pre_or_post choice\n"));
      fun2 = NULL;
      sgs_nums = NULL;
   }
   return(status);
}

/* ------------------------------------------------------------------------- */
/* generate some other Gauss Seidel smoother (Ray, what is it ?)             */
/* ------------------------------------------------------------------------- */

int ML_Gen_SmootherGSextra( ML *ml , int nl, int pre_or_post, int ntimes, 
                            double omega, int Nextra, int extra[])
{
   int (*fun)(void *, int, double *, int, double *);
   int start_level, end_level, i, status, Nrows, count;
   double *nums, **sgs_nums = NULL, *num2;
   ML_Operator *Amat;
   struct ML_CSR_MSRdata *ptr = NULL;
   double *val = NULL, temp_omega;
   int *bindx;
   void (*fun2)(void *) = NULL;


   if (nl == ML_ALL_LEVELS) { start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_SmootherGSextra: cannot set smoother on level %d\n",start_level);
      return 1;
   }
	
   if (omega == ML_DEFAULT) omega = 1.;
   Amat = &(ml->Amat[nl]);
   fun  = ML_Smoother_SGS;

   if (Amat->getrow->external == MSR_getrows){
      ptr   = (struct ML_CSR_MSRdata *) Amat->data;
      val   = ptr->values;
      bindx = ptr->columns;
   }
#ifdef AZTEC
   else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif
   if (val != NULL) {
      fun = ML_MSR_SGSextra;
      sgs_nums = (double **) ML_allocate( sizeof(double)*4);
      Nrows    = Amat->getrow->Nrows;
      nums     = (double *) ML_allocate( Nrows * sizeof(double));
      num2     = (double *) ML_allocate( Nrows * sizeof(double));
      sgs_nums[0] = nums;
      sgs_nums[1] = num2;
sgs_nums[2] = (double *) ML_allocate( 1 * sizeof(double) );
sgs_nums[2][0] = (double) Nextra;
sgs_nums[3] = (double *) extra;
      for (i = 0; i < Nrows; i++) { 
         count = 0;
/*
         for (j = bindx[i]; j < bindx[i+1]; j++) 
            if (bindx[j] >= Nrows) count++;
*/

         if (bindx[i] != bindx[i+1]) 
            temp_omega = omega*(1.0 - .5*((double) count)/ ((double) (bindx[i+1]-bindx[i])));
         else temp_omega = 1.;

         num2[i] = 1. - temp_omega;
         if (val[i] != 0.0) nums[i] = temp_omega/val[i];
         else { nums[i] = 0.0; num2[i] = 1.; }
      }
      fun2 = ML_MSR_GSextra_Clean;
   }

	
   if (pre_or_post == ML_PRESMOOTH) 
      for (i = start_level; i <= end_level; i++) {
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL, sgs_nums,
				  fun, NULL, ntimes, omega);
	 ml->pre_smoother[i].data_destroy = fun2;
   }
   else if (pre_or_post == ML_POSTSMOOTH)
      for (i = start_level; i <= end_level; i++) {
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL, sgs_nums,
                             fun, NULL, ntimes, omega);
	 ml->post_smoother[i].data_destroy = fun2;
      }
   else return(pr_error("Print unknown pre_or_post choice\n"));
   return(status);
}

/* ------------------------------------------------------------------------- */
/* generate the ordered symmetric Gauss Seidel smoother                      */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_OrderedSymGaussSeidel(ML *ml , int nl, int pre_or_post,
		 		          int ntimes, double omega)
{
   int (*fun)(void *, int, double *, int, double *);
   int start_level, end_level, i, status, *ordering;
#ifdef ML_TIMING
   double t0;
#endif

   if (nl == ML_ALL_LEVELS) { start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_Smoother_OrderedSymGaussSeidel: cannot set smoother on level %d\n",start_level);
      return 1;
   }
	
   fun = ML_Smoother_OrderedSGS;
   if (omega == ML_DEFAULT) omega = 1.;
	
   if (pre_or_post == ML_PRESMOOTH) 
      for (i = start_level; i <= end_level; i++) {
#ifdef ML_TIMING
         t0 = GetClock();
#endif
         ML_Smoother_Gen_Ordering(&(ml->Amat[i]), &ordering);
         ml->pre_smoother[i].data_destroy = ML_Smoother_Clean_OrderedSGS;
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL, 
                         (void *) ordering, fun, NULL, ntimes, omega);
#ifdef ML_TIMING
         ml->pre_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->pre_smoother[i].build_time;
#endif
      }
   else if (pre_or_post == ML_POSTSMOOTH)
      for (i = start_level; i <= end_level; i++) {
#ifdef ML_TIMING
         t0 = GetClock();
#endif
         ML_Smoother_Gen_Ordering(&(ml->Amat[i]), &ordering);
         ml->post_smoother[i].data_destroy = ML_Smoother_Clean_OrderedSGS;
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL, 
                         (void *) ordering, fun, NULL, ntimes, omega);
#ifdef ML_TIMING
         ml->post_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->post_smoother[i].build_time;
#endif
      }
   else return(pr_error("Print unknown pre_or_post choice\n"));
   return(status);
}

/* ------------------------------------------------------------------------- */
/* generate the block Gauss Seidel smoother (fixed size block)               */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_BlockGaussSeidel(ML *ml , int nl, int pre_or_post,
                                     int ntimes, double omega, int blocksize)
{
   int            (*fun)(void *, int, double *, int, double *);
   ML_Sm_BGS_Data *data;
   int            start_level, end_level, i, status;
#ifdef ML_TIMING
   double         t0;
   t0 = GetClock();
#endif

   if (nl == ML_ALL_LEVELS) { start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_Smoother_BlockGaussSeidel: cannot set smoother on level %d\n",start_level);
      return 1;
   }
	
   fun = ML_Smoother_BlockGS;
   if (omega == ML_DEFAULT) omega = 1.;

   if (pre_or_post == ML_PRESMOOTH) {
      for (i = start_level; i <= end_level; i++) {
         ML_Smoother_Create_BGS_Data(&data);
	 ML_Smoother_Gen_BGSFacts(&data, &(ml->Amat[i]), blocksize);
	 ml->pre_smoother[i].data_destroy = ML_Smoother_Clean_BGS_Data;
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL,
		                 (void *) data, fun, NULL, ntimes, omega);
#ifdef ML_TIMING
         ml->pre_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->pre_smoother[i].build_time;
#endif
      }
   }
   else if (pre_or_post == ML_POSTSMOOTH) {
      for (i = start_level; i <= end_level; i++) {
         ML_Smoother_Create_BGS_Data(&data);
	 ML_Smoother_Gen_BGSFacts(&data, &(ml->Amat[i]), blocksize);
	 ml->post_smoother[i].data_destroy = ML_Smoother_Clean_BGS_Data;
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL,
			      (void *) data, fun, NULL, ntimes, omega);
#ifdef ML_TIMING
         ml->post_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->post_smoother[i].build_time;
#endif
      }
   }
   else if (pre_or_post == ML_BOTH) {
      for (i = start_level; i <= end_level; i++) {
         ML_Smoother_Create_BGS_Data(&data);
	 ML_Smoother_Gen_BGSFacts(&data, &(ml->Amat[i]), blocksize);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL,
		                 (void *) data, fun, NULL, ntimes, omega);
#ifdef ML_TIMING
         ml->pre_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->pre_smoother[i].build_time;
#endif
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL,
			      (void *) data, fun, NULL, ntimes, omega);
      }
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
   return(status);
}

/* ------------------------------------------------------------------------- */
/* generate the variable block Jacobi smoother                               */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockJacobi( ML *ml , int nl, int pre_or_post, 
                 int ntimes, double omega, int nBlocks, int *blockIndices)
{
   int            (*fun)(void *, int, double *, int, double *);
   double         myomega;
   ML_Sm_BGS_Data *data;
	
   if (nl == ML_ALL_LEVELS) { 
      printf("ML_Gen_Smoother_VBlockJacobi: ML_ALL_LEVELS not allowed\n");
      return 1;
   }
   if (nl < 0) {
      printf("ML_Gen_Smoother_VBlockJacobi: cannot set smoother on level %d\n",nl);
      return 1;
   }
   fun = ML_Smoother_VBlockJacobi;
   if (omega == ML_DEFAULT) myomega = .5;
   else                     myomega = omega;
	
   ML_Smoother_Create_BGS_Data(&data);

   ML_Smoother_Gen_VBGSFacts(&data, &(ml->Amat[nl]), nBlocks, blockIndices);

   if (pre_or_post == ML_PRESMOOTH) 
         return(ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                        (void *) data, fun, NULL, ntimes, myomega));
   else if (pre_or_post == ML_POSTSMOOTH)
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, myomega));
   else if (pre_or_post == ML_BOTH) {
      ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                        (void *) data, fun, NULL, ntimes, myomega);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, myomega));
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the variable block Gauss Seidel smoother                         */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockSymGaussSeidel( ML *ml , int nl, int pre_or_post, 
                      int ntimes, double omega, int nBlocks, int *blockIndices)
{
   int            (*fun)(void *, int, double *, int, double *);
   ML_Sm_BGS_Data *data;
	
   if (nl == ML_ALL_LEVELS) { 
      printf("ML_Gen_Smoother_VBlockSymGaussSeidel: ML_ALL_LEVELS not allowed\n");
      return 1;
   }
   if (nl < 0) {
      printf("ML_Gen_Smoother_VBlockSymGaussSeidel: cannot set smoother on level %d\n",nl);
      return 1;
   }
   fun = ML_Smoother_VBlockSGS;
	
   ML_Smoother_Create_BGS_Data(&data);

   ML_Smoother_Gen_VBGSFacts(&data, &(ml->Amat[nl]), nBlocks, blockIndices);

   if (pre_or_post == ML_PRESMOOTH) 
         return(ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                        (void *) data, fun, NULL, ntimes, omega));
   else if (pre_or_post == ML_POSTSMOOTH)
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, omega));
   else if (pre_or_post == ML_BOTH) {
      ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                        (void *) data, fun, NULL, ntimes, omega);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, omega));
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the variable block Gauss Seidel smoother (sequential)            */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockSymGaussSeidelSequential( ML *ml , int nl, 
     int pre_or_post, int ntimes, double omega, int nBlocks, int *blockIndices)
{
   int            (*fun)(void *, int, double *, int, double *);
   ML_Sm_BGS_Data *data;
	
   if (nl == ML_ALL_LEVELS) { 
      printf("ML_Gen_Smoother_VBlockSymGaussSeidelSequential: ML_ALL_LEVELS not allowed\n");
      return 1;
   }
   if (nl < 0) {
      printf("ML_Gen_Smoother_VBlockSymGaussSeidelSequential: cannot set smoother on level %d\n",nl);
      return 1;
   }
   fun = ML_Smoother_VBlockSGSSequential;
	
   ML_Smoother_Create_BGS_Data(&data);

   ML_Smoother_Gen_VBGSFacts(&data, &(ml->Amat[nl]), nBlocks, blockIndices);

   if (pre_or_post == ML_PRESMOOTH) 
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, omega));
   else if (pre_or_post == ML_POSTSMOOTH)
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, omega));
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the variable block Jacobi smoother with Krylov(not debugged yet) */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockKrylovJacobi( ML *ml , int nl, int pre_or_post, 
                 int ntimes, double omega, int nBlocks, int *blockIndices)
{
   int            (*fun)(void *, int, double *, int, double *);
   double         myomega;
   ML_Sm_BGS_Data *data;
	
   if (nl == ML_ALL_LEVELS) { 
      printf("ML_Gen_Smoother_VBlockKrylovJacobi: ML_ALL_LEVELS not allowed\n");
      return 1;
   }
   if (nl < 0) {
      printf("ML_Gen_Smoother_VBlockKrylovJacobi: cannot set smoother on level %d\n",nl);
      return 1;
   }
   fun = ML_Smoother_VBlockKrylovJacobi;
   if (omega == ML_DEFAULT) myomega = .5;
   else                     myomega = omega;
	
   ML_Smoother_Create_BGS_Data(&data);

   ML_Smoother_Gen_VBGSFacts(&data, &(ml->Amat[nl]), nBlocks, blockIndices);

   if (pre_or_post == ML_PRESMOOTH) 
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, myomega));
   else if (pre_or_post == ML_POSTSMOOTH)
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, myomega));
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the overlapped domain decomposition with ILUT smoother           */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_OverlappedDDILUT( ML *ml , int nl, int pre_or_post )
{
   int             (*fun)(void *, int, double *, int, double *);
   int             total_recv_leng, *recv_lengths, *int_buf, *map, *map2; 
   int             offset;
   double          *dble_buf;
   ML_Sm_ILUT_Data *data;
   ML_Operator     *Amat;
   ML_Comm         *comm;
	
   /* ---------------------------------------------------------------- */
   /* initialize the ILUT data object                                  */
   /* ---------------------------------------------------------------- */
   if (nl == ML_ALL_LEVELS) { 
      printf("ML_Gen_Smoother_OverlappedDDILUT: ML_ALL_LEVELS not allowed\n");
      return 1;
   }
   if (nl < 0) {
      printf("ML_Gen_Smoother_OverlappedDDILUT: cannot set smoother on level %d\n",nl);
      return 1;
   }

   fun = ML_Smoother_OverlappedILUT;
	
   comm = ml->comm;
   Amat = &(ml->Amat[nl]);
   ML_Smoother_Create_ILUT_Data( &data );
   data->fillin    = 1;
   data->threshold = 1.0e-8;

   /* ---------------------------------------------------------------- */
   /* send the lengths of each row to remote processor at the end,     */
   /* additional row information should be given in total_recv_leng,   */
   /* recv_lengths, int_buf, dble_buf                                  */
   /* ---------------------------------------------------------------- */

   ML_Smoother_ComposeOverlappedMatrix(Amat, comm, &total_recv_leng, 
              &recv_lengths, &int_buf, &dble_buf, &map, &map2, &offset);

   /* ---------------------------------------------------------------- */
   /* use the local matrix row and the off-processor rows to compose   */
   /* ILUT preconditioner                                              */
   /* ---------------------------------------------------------------- */

   ML_Smoother_ILUTDecomposition(data,Amat,comm,total_recv_leng,recv_lengths, 
                                 int_buf, dble_buf, map, map2, offset);
   if ( map  != NULL ) free(map);
   if ( map2 != NULL ) free(map2);
   if ( int_buf != NULL ) free(int_buf);
   if ( dble_buf != NULL ) free(dble_buf);
   if ( recv_lengths != NULL ) free(recv_lengths);

   /* ---------------------------------------------------------------- */
   /* set it up as smoother                                            */
   /* ---------------------------------------------------------------- */

   if (pre_or_post == ML_PRESMOOTH) 
         return(ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                        (void *) data, fun, NULL, 1, 0.0));
   else if (pre_or_post == ML_POSTSMOOTH)
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, 1, 0.0));
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the variable block additive Schwarz smoother                     */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockAdditiveSchwarz(ML *ml , int nl, int pre_or_post,
                                          int ntimes, int length, int *blkinfo)
{
   int                (*fun)(void *, int, double *, int, double *);
   int                total_recv_leng, *recv_lengths, *int_buf, *map, *map2; 
   int                i, maxblk, offset;
   double             *dble_buf;
   ML_Sm_Schwarz_Data *data;
   ML_Operator        *Amat;
   ML_Comm            *comm;
	
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
   data->blk_info = (int *) malloc(data->Nrows * sizeof(int));
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
   if ( map  != NULL ) free(map);
   if ( map2 != NULL ) free(map2);
   if ( int_buf != NULL ) free(int_buf);
   if ( dble_buf != NULL ) free(dble_buf);
   if ( recv_lengths != NULL ) free(recv_lengths);

   /* ---------------------------------------------------------------- */
   /* set it up as smoother                                            */
   /* ---------------------------------------------------------------- */

   if (pre_or_post == ML_PRESMOOTH) 
         return(ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                        (void *) data, fun, NULL, ntimes, 0.0));
   else if (pre_or_post == ML_POSTSMOOTH)
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, 0.0));
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the variable block additive Schwarz smoother                     */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockMultiplicativeSchwarz(ML *ml , int nl, int pre_or_post,
                                         int ntimes, int length, int *blkinfo )
{
   int                (*fun)(void *, int, double *, int, double *);
   int                total_recv_leng, *recv_lengths, *int_buf, *map, *map2; 
   int                i, maxblk, offset;
   double             *dble_buf;
   ML_Sm_Schwarz_Data *data;
   ML_Operator        *Amat;
   ML_Comm            *comm;
	
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
   data->blk_info = (int *) malloc(data->Nrows * sizeof(int));
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
   if ( map  != NULL ) free(map);
   if ( map2 != NULL ) free(map2);
   if ( int_buf != NULL ) free(int_buf);
   if ( dble_buf != NULL ) free(dble_buf);
   if ( recv_lengths != NULL ) free(recv_lengths);

   /* ---------------------------------------------------------------- */
   /* set it up as smoother                                            */
   /* ---------------------------------------------------------------- */

   if (pre_or_post == ML_PRESMOOTH) 
         return(ML_Smoother_Set(&(ml->pre_smoother[nl]), ML_INTERNAL, 
                        (void *) data, fun, NULL, ntimes, 0.0));
   else if (pre_or_post == ML_POSTSMOOTH)
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), ML_INTERNAL, 
                             (void *) data, fun, NULL, ntimes, 0.0));
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the sparse approximate inverse smoother */
/* ------------------------------------------------------------------------- */

#ifdef PARASAILS
#include "Matrix.h"
#include "ParaSails.h"

#endif

int ML_Gen_Smoother_ParaSails(ML *ml, int nl, int pre_or_post, int ntimes,
   int sym, double thresh, int num_levels, double filter, int parasails_loadbal,
   int parasails_factorized)
{
#ifdef PARASAILS
   int            (*fun1)(void *, int, double *, int, double *);
   int            (*fun2)(void *, int, double *, int, double *);
   int            start_level, end_level, i, status;
   int            row, start_row, end_row, row_length;

   Matrix *mat;
   struct widget { int parasails_factorized; ParaSails *ps;} *widget;
   ParaSails *ps;
   int j;

#ifdef ML_TIMING
   double         t0;
   t0 = GetClock();
#endif


   if (nl == ML_ALL_LEVELS) { start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_Smoother_ParaSails: cannot set smoother on level %d\n",start_level);
      return 1;
   }
	
   if (sym)
      fun1 = fun2 = ML_Smoother_ParaSailsSym;
   else {
      fun1 = ML_Smoother_ParaSails;
      fun2 = ML_Smoother_ParaSailsTrans;
   }

      for (i = start_level; i <= end_level; i++) {

	 int nrows = ml->Amat[i].outvec_leng;
	 int local_row, allocated_space, *ml_indices;
         double *ml_values;
	 start_row = ML_gpartialsum_int(nrows, ml->comm);
	 end_row   = start_row + nrows - 1; 

	 mat = MatrixCreate(ml->comm->USR_comm, start_row, end_row);

	 ML_create_unique_id(ml->Amat[i].invec_leng,
	                     &(ml->Amat[i].getrow->loc_glob_map),
	                     ml->Amat[i].getrow->pre_comm, ml->comm);
	 ml->Amat[i].getrow->use_loc_glob_map = ML_YES;
	 allocated_space = 10; 
	 ml_indices = (int    *) ML_allocate(sizeof(int)*allocated_space);
	 ml_values  = (double *) ML_allocate(sizeof(double)*allocated_space);

	 for (row=start_row; row<=end_row; row++) {
	    /* global indices */

            local_row = row - start_row;
            ML_get_matrix_row(&(ml->Amat[i]), 1, &local_row,
                                &allocated_space, &ml_indices, 
				&ml_values, &row_length, 0);
/*
for (j = 0; j < row_length; j++)
   printf("A(%d,%d) = %e;\n",row+1,ml_indices[j]+1,ml_values[j]);
*/

	    MatrixSetRow(mat, row, row_length, ml_indices, ml_values);
	 }
	 ML_free(ml->Amat[i].getrow->loc_glob_map); 
         ml->Amat[i].getrow->loc_glob_map = NULL;
	 ml->Amat[i].getrow->use_loc_glob_map = ML_NO;

	 MatrixComplete(mat);

	 /* nonsymmetric preconditioner */

         widget = (struct widget *) ML_allocate(sizeof(struct widget));
	 ps  = ParaSailsCreate(ml->comm->USR_comm, start_row, end_row, 
	    parasails_factorized);
         ps->loadbal_beta = parasails_loadbal;
	 ParaSailsSetupPattern(ps, mat, thresh, num_levels);
	 ParaSailsStatsPattern(ps, mat);
	 ParaSailsSetupValues(ps, mat, filter);
	 ParaSailsStatsValues(ps, mat);

	 /* we can destroy the matrix now */
	 MatrixDestroy(mat);

	 /* ml->post_smoother[i].data_destroy = ML_Smoother_Clean_ParaSails; */

/* Turn on temporarily */
         widget->parasails_factorized = parasails_factorized;
         widget->ps                   = ps;
         
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), ML_INTERNAL,
			      (void *) widget, fun1, NULL, ntimes, 0.0);

	 status = ML_Smoother_Set(&(ml->pre_smoother[i]), ML_INTERNAL,
			      (void *) widget, fun2, NULL, ntimes, 0.0);

#ifdef ML_TIMING
         ml->post_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->post_smoother[i].build_time;
#endif
      }

   /* note: in free, post and pre are the same */

   return(status);
#else
   printf("ParaSails not linked\n");
   return(1);
#endif
}

/* ************************************************************************* */
/* functions to generate Galerkin coarse grid operator                       */
/* ------------------------------------------------------------------------- */

int ML_Gen_AmatrixRAP(ML *ml, int parent_level, int child_level)
{
   ML_Operator *Amat, *Rmat, *Pmat;
   int i, output_level;

#ifdef ML_TIMING
   double t0;
   t0 = GetClock();
#endif

   output_level = ml->output_level;
   i   = parent_level;
   Amat = &(ml->Amat[parent_level]);
   Rmat = &(ml->Rmat[parent_level]);
   Pmat = &(ml->Pmat[child_level]);

   if (Amat->matvec->ML_id == ML_EMPTY) {
      if (output_level > 3)
      printf("Warning: No Amat matvec on grid %d (where finest = 0).\n\
		can not check Amat's getrow\n",i);
   }
/*
   else ML_Operator_Check_Getrow( Amat,ml->comm,i,"Amat");
*/

   if (Amat->getrow->ML_id == ML_EMPTY)
      pr_error("Error: No A matrix getrow on grid %d : \
                       can not do ML_Gen_Amatrix_RAP.\n",i);

   if ((Amat->getrow->pre_comm  == NULL) && 
       (Amat->getrow->post_comm == NULL) && (ml->comm->ML_nprocs > 1) ) {
       if ((ml->comm->ML_mypid == 0) && (output_level > 3)) {
         printf("Warning:No communication information given with Amat's \n");
         printf("\tgetrow on level %d (finest = 0).!!!!\n",i);
       }
   }

   if (Rmat->matvec->ML_id == ML_EMPTY) {
      if (output_level > 3)
      printf("Warning: No Rmat matvec on grid %d (where finest = 0).\n\
		can not check Rmat's getrow\n",i);
   }
/*
   else ML_Operator_Check_Getrow( Rmat,ml->comm,i,"Rmat");
*/

   if (Rmat->getrow->ML_id == ML_EMPTY)
      pr_error("Error: No R matrix getrow on grid %d : \n\
                       can not do ML_Gen_AmatrixRAP.\n",i);

   if ((Rmat->getrow->pre_comm  == NULL) && 
       (Rmat->getrow->post_comm == NULL) && (ml->comm->ML_nprocs > 1) ) {
       if ((ml->comm->ML_mypid == 0) && (output_level > 3)) {
         printf("Warning:No communication information given with Rmat's \n");
         printf("\tgetrow on level %d (finest = 0).!!!!\n",i);
       }
   }

   if (Pmat->matvec->ML_id == ML_EMPTY) {
      if (output_level > 3) 
      printf("Warning: No Pmat matvec on grid %d (where finest = 0).\n\
		can not check Pmat's getrow\n",i);
   }
/*
   else ML_Operator_Check_Getrow(Pmat,ml->comm,i,"Pmat");
*/

   if (Pmat->getrow->ML_id == ML_EMPTY)
      pr_error("Error: No P matrix getrow on grid %d : \n\
                       can not do ML_Gen_AmatrixRAP.\n",i);

   if ((Pmat->getrow->pre_comm  == NULL) && 
       (Pmat->getrow->post_comm == NULL) && (ml->comm->ML_nprocs > 1) ) {
       if ((ml->comm->ML_mypid == 0) && (output_level > 3)) {
         printf("Warning:No communication information given with Pmat's \n");
         printf("\tgetrow on level %d (finest = 0).!!!!\n",i);
       }
   }

   ML_rap(ml, &(ml->Rmat[parent_level]), &(ml->Amat[parent_level]), 
          &(ml->Pmat[child_level]), &(ml->Amat[child_level]), 
	  child_level, ml->Pmat[child_level].invec_leng);
#ifdef ML_TIMING
   ml->Amat[child_level].build_time = GetClock() - t0;
   ml->timing->total_build_time   += ml->Amat[child_level].build_time;
#endif
   return(1);
}

/* ************************************************************************* */
/* functions to set up boundary                                              */
/* ------------------------------------------------------------------------- */

int ML_Set_BoundaryTypes(ML *ml, int level, int type, int n, int *data)
{
   ML_BdryPts *ml_bc;

   if ( type != ML_BDRY_DIRICHLET )
   {
      printf("ML_Set_BoundaryTypes error : type not supported.\n");
      exit(1);
   }
   ml_bc = ml->SingleLevel[level].BCs;
   ML_BdryPts_Load_Dirichlet_Grid(ml_bc, n, data);
   return 0;
}

/* ************************************************************************* */
/* functions to set the mapper                                               */
/* ------------------------------------------------------------------------- */

int ML_Set_EqnToGridMapFunc(ML *ml, int level, int fleng, int tleng, 
                           void *data, int (*func)(void*,double*,double*) )
{
   ML_Mapper_SetFunc( ml->SingleLevel[level].eqn2grid,fleng,tleng,func);
   ML_Mapper_SetData( ml->SingleLevel[level].eqn2grid, 
                      (void*)ml->SingleLevel[level].Grid->Grid);
   return 0;
}

int ML_Set_GridToEqnMapFunc(ML *ml, int level, int fleng, int tleng, 
                           void *data, int (*func)(void*,double*,double*) )
{
   ML_Mapper_SetFunc( ml->SingleLevel[level].grid2eqn,fleng,tleng,func);
   ML_Mapper_SetData( ml->SingleLevel[level].grid2eqn,
                      (void*)ml->SingleLevel[level].Grid->Grid);
   return 0;
}

/* ************************************************************************* */
/* Setting up a given solver                                                 */
/* ------------------------------------------------------------------------- */

int ML_Gen_Solver(ML *ml, int scheme, int finest_level, int coarsest_level)
{
   int        i, j, level, leng1, leng2, leng3, *itmp3;
   int        output_level, t1, t2;
   double     *dtmp1, *dtmp2;
   ML_1Level  *current_level, *temp;

   ml->ML_scheme         = scheme;
   ml->ML_finest_level   = finest_level;
   ml->ML_coarsest_level = coarsest_level;
   output_level          = ml->output_level;

   if (output_level > 1) {
      if (ml->comm->USR_sendbytes  == NULL) {
         printf("Warning: Machine's communication platform (e.g. MPI) not\n");
         printf("         set. Assume we are running in serial.\n");
      }
   }


   t1                    = finest_level - coarsest_level;
   if ( t1 < 0 ) t1 = - t1;

   /*
   if ( (t1+1) != nlevels) {
      pr_error("ML_Make_Solver Error : inconsistent level information.\n");
      exit(1);
   }
   */
      
   current_level = &(ml->SingleLevel[finest_level]);
   level = finest_level;
   i = 0;
   while (current_level != NULL) {
      if (current_level->Amat->matvec->ML_id == ML_EMPTY &&
          level != coarsest_level) {
         pr_error("Error: No A matrix on grid %d.\n",level);
      }

      if ((current_level->Amat->getrow->pre_comm  == NULL) && 
          (current_level->Amat->getrow->post_comm == NULL) && 
          (current_level->Amat->getrow->ML_id != ML_EMPTY) &&
          (ml->comm->ML_nprocs > 1) ) {
         if (ml->comm->ML_mypid == 0) {
         printf("Warning:No communication information given with Amat's \n");
         printf("\tgetrow on level %d (finest = 0).!!!!\n",level); }
      }
      else {
         /* This does not work if Nrows in a processor = 0 */
/*
            ML_Operator_Check_Getrow(current_level->Amat,ml->comm,level,"Amat");
*/
      }

      temp = current_level->Rmat->to;

      if (temp != NULL) {
         if (current_level->Rmat->matvec->ML_id == ML_EMPTY)
            pr_error("Error: No R matvec on grid %d.\n",level);
/*
         ML_Operator_Check_Getrow(current_level->Rmat,ml->comm,level,"Rmat");
*/
         if (level != finest_level && 
             current_level->Pmat->matvec->ML_id == ML_EMPTY)
            pr_error("Error: No P matvec on grid %d.\n",level);
/*
         ML_Operator_Check_Getrow(current_level->Pmat,ml->comm,level, "Pmat");
*/
      }

      if ( (current_level->pre_smoother->smoother->ML_id == ML_INTERNAL) &&
           (current_level->pre_smoother->smoother->internal==ML_Smoother_Jacobi)) 
      {
         if ((temp == NULL) && (current_level->csolve->func->ML_id == ML_EMPTY))
         {
            if (current_level->pre_smoother->ntimes == ML_NOTSET) {
               current_level->pre_smoother->ntimes = ML_CONVERGE;
               current_level->pre_smoother->tol    = 1.0e-10;
               if ((output_level > 3) && (ml->comm->ML_mypid == 0))  {
                  printf("Iterating Jacobi on grid %d until\n",level);
                  printf("convergence.  This could be very very slow.\n");
               }
            }
         }
/*
         else if ((ml->comm->ML_mypid == 0) && (output_level > 3)) {
            printf("Warning: Using the ML Jacobi smoother on grid %d\n",level);
            printf("         This can be inefficient.        \n");
         }
*/
      }
      if (current_level->pre_smoother->ntimes == ML_NOTSET) 
            current_level->pre_smoother->ntimes = 2;
      if (temp != NULL) {
         t1 = current_level->Amat->outvec_leng;
         t2 = temp->Amat->outvec_leng;
         ML_gsum_vec_int(&t1, &j, 1, ml->comm);
         ML_gsum_vec_int(&t2, &j, 1, ml->comm);
         if (t2 >= t1) {
            if (ml->comm->ML_mypid == 0) 
               pr_error("Error: Grid %d (where finest = 0) has %d unknowns \
                       and restricts to a grid with %d (i.e. more) unknowns.\n",
	                i, t1, t2);
            else pr_error("");
         }
      }
      i++;
      if ( coarsest_level < finest_level ) level--; else level++;

      if ( ML_BdryPts_Check_Dirichlet_Grid(current_level->BCs) == 1 &&
           ML_Mapper_Check(current_level->grid2eqn) == 1 )
      {
         ML_Mapper_GetLength(current_level->grid2eqn, &leng1, &leng2);
         dtmp1 = (double*) malloc( leng1 * sizeof(double) );
         dtmp2 = (double*) malloc( leng2 * sizeof(double) );
         ML_BdryPts_Get_Dirichlet_Grid_Info(current_level->BCs,&leng3,&itmp3);
         for ( j = 0; j < leng1; j++ ) dtmp1[j] = 0.0;
         for ( j = 0; j < leng2; j++ ) dtmp2[j] = 0.0;
         for ( j = 0; j < leng3; j++ ) dtmp1[itmp3[j]] = 1.0;
         ML_Mapper_Apply(current_level->grid2eqn, dtmp1, dtmp2);
         leng1 = 0;
         for ( j = 0; j < leng2; j++ ) if ( dtmp2[j] == 1.0 ) leng1++;
         itmp3 = (int*) malloc( leng1 * sizeof(int) );
         leng1 = 0;
         for ( j = 0; j < leng2; j++ ) 
            if ( dtmp2[j] == 1.0 ) itmp3[leng1++] = j;
         ML_BdryPts_Load_Dirichlet_Eqn(current_level->BCs, leng1, itmp3);
         free( itmp3 );
      } else {
         ML_BdryPts_Copy_Dirichlet_GridToEqn(current_level->BCs);
      }    
      current_level = temp;
   }
   if ((ml->comm->ML_mypid == 0) && (output_level > 0)) 
      printf("Total number of levels = %d\n",i);
  
   if ((output_level > 5) && (ml->comm->ML_mypid == 0)) {
      if (i == 1) printf("Warning: Only a one level multilevel scheme!!\n");
      /*
      if (i > 2) 
         printf("Warning: Its best to start with 2 levels.\n");
      */
   }
   
   if ( finest_level > coarsest_level ) {
      for ( i = coarsest_level; i < finest_level; i++ )
         ml->Pmat[i].bc = &(ml->BCs[i+1]);
   } else {
      for ( i = coarsest_level; i > finest_level; i-- )
         ml->Pmat[i].bc = &(ml->BCs[i-1]);
   }
   return 0;
}

/* ************************************************************************* */
/* ML_Iterate : call iteration function                                      */
/*-------------------------------------------------------------------------- */

int ML_Iterate(ML *ml, double *sol, double *rhs)
{
   int  i = 0, count = 0;
   double res_norm, prev_res_norm = -1.0, reduction, r0/*, old_reduction = 10.*/;

   reduction = 1.;

   if ((ml->res_output_freq > 0) && (ml->comm->ML_mypid == 0) ) {
      printf("Iter       ||res_i||_2**    ||res_i||/||res_i+1||\n");
      count = ml->res_output_freq - 1;
   }
   else count = 0;

   while ( reduction >= ml->tolerance && i < ml->max_iterations ) 
   {
      res_norm = ML_Cycle_MGV( &(ml->SingleLevel[ml->ML_finest_level]), sol, 
				rhs, ML_NONZERO,ml->comm, ML_COMPUTE_RES_NORM);
      count++;
      i++;
      if (count == ml->res_output_freq) {
         if ((ml->comm->ML_mypid == 0) && (i == 1))
            printf("%4d       %10.3e\n",i,res_norm);
         if ((ml->comm->ML_mypid == 0) && (prev_res_norm != -1.))
            printf("%4d       %10.3e           %10.3e\n",i,res_norm,
		   res_norm/prev_res_norm);
         count = 0;
      }
      if (i == 1) r0 = res_norm + 1.e-15;
      prev_res_norm = res_norm + 1.e-15;
      reduction = res_norm/r0;
/*
      if (reduction >= old_reduction) {
         if (ml->comm->ML_mypid == 0) printf("**** Stagnation  ******\n");
         reduction = -1.;
      }
      if ( (i%5) == 1) old_reduction = reduction;
*/
   }

   if ((ml->res_output_freq > 0) && (ml->comm->ML_mypid == 0) )
      printf("\n**Residual norm taken after multigrid pre_smoothing step.\n\n");
   return 0;
}

/*****************************************************************************/
/* solve using V-cycle multigrid                                             */
/*-------------------------------------------------------------------------- */

int ML_Solve_MGV( ML *ml , double *din, double *dout)
{
   int    i, leng, dir_leng, *dir_list, k, level;
   double *diag, *scales, *din_temp;

   /* ------------------------------------------------------------ */
   /* initially set the solution to be all 0           	           */
   /* ------------------------------------------------------------ */

   level = ml->ML_finest_level;
   leng = ml->Amat[level].outvec_leng;
   for ( i = 0; i < leng; i++ ) dout[i] = 0.0;
   din_temp = (double*) malloc( leng * sizeof(double) );

   /* ------------------------------------------------------------ */
   /* on the fixed boundaries, set them to be each to the solution */
   /* ------------------------------------------------------------ */

   ML_BdryPts_Get_Dirichlet_Eqn_Info(&(ml->BCs[level]),&dir_leng,&dir_list);
   if ( dir_leng != 0 )
   {
      if (ml->Amat[level].diagonal != NULL) { 
         ML_DVector_GetDataPtr(ml->Amat[level].diagonal,&diag);
         for ( i = 0; i < dir_leng; i++ ) {
            k = dir_list[i]; 
            dout[k] = din[k] / diag[k];
         } 
      } else {
         diag = NULL;
         for ( i = 0; i < dir_leng; i++ ) {
            k = dir_list[i]; 
            dout[k] = din[k];
         }
      }
   }

   /* ------------------------------------------------------------ */
   /* normalization                                                */
   /* ------------------------------------------------------------ */

   ML_DVector_GetDataPtr(&(ml->Amat_Normalization[level]), &scales) ;

/* watch out for this !!!!! */
scales = NULL;

   if ( scales != NULL ) {
      for ( i = 0; i < leng; i++ ) din_temp[i] = din[i] / scales[i];
   } else {
      scales = NULL;
      for ( i = 0; i < leng; i++ ) din_temp[i] = din[i];
   }

   /* ------------------------------------------------------------ */
   /* call MG v-cycle                                              */
   /* ------------------------------------------------------------ */

   ML_Cycle_MGV(&(ml->SingleLevel[ml->ML_finest_level]), dout, din_temp, 
                ML_ZERO, ml->comm, ML_NO_RES_NORM);

   free(din_temp);
   return 0;
}

void ML_Solve_SmootherDestroy(void *data) 
{
   ML_Destroy( (ML **) &data);
}
int ML_Solve_Smoother(void *data, int isize, double *x, int osize, double *rhs)
{
   ML *ml;
   int n, i;
   double *res,*tmp;

   ml = (ML *) data;
   n = ml->Amat[0].invec_leng;
   tmp  = (double *) ML_allocate(n*sizeof(double));
   res  = (double *) ML_allocate(n*sizeof(double));
   if (res == NULL) pr_error("swillie: out of space\n");

   ML_Operator_Apply(&(ml->Amat[0]), n, x, n, res);
   for (i = 0; i < n; i++) res[i] = rhs[i] - res[i];
   for (i = 0; i < n; i++) tmp[i] = 0.;

   ML_Solve_MGV( ml, res, tmp);

   for (i = 0; i < n; i++) x[i] += tmp[i];
   ML_free(res);
   ML_free(tmp);
   return 0;
}

/*****************************************************************************/
/* segregated solve                                                          */
/*-------------------------------------------------------------------------- */

int ML_Seg_Solve( ML *ml , double *din, double *dout)
{
   int    i, leng, dir_leng, *dir_list, k, level;
   double *diag, *scales, *din_temp;

   /* ------------------------------------------------------------ */
   /* initially set the right hand side to be all 0                */
   /* ------------------------------------------------------------ */

   level = ml->ML_finest_level;
   leng = ml->Amat[level].outvec_leng;
   for ( i = 0; i < leng; i++ ) dout[i] = 0.0;
   din_temp = (double*) malloc( leng * sizeof(double) );

   /* ------------------------------------------------------------ */
   /* on the fixed boundaries, set them to be each to the solution */
   /* ------------------------------------------------------------ */

   ML_BdryPts_Get_Dirichlet_Eqn_Info(&(ml->BCs[level]),&dir_leng,&dir_list);
   if ( dir_leng != 0 )
   {
      if (ml->Amat[level].diagonal != NULL) { 
         ML_DVector_GetDataPtr(ml->Amat[level].diagonal,&diag);
         for ( i = 0; i < dir_leng; i++ ) {
            k = dir_list[i]; 
            dout[k] = din[k] / diag[k];
         } 
      } else {
         diag = NULL;
         for ( i = 0; i < dir_leng; i++ ) {
            k = dir_list[i]; 
            dout[k] = din[k];
         }
      }
   }

   /* ------------------------------------------------------------ */
   /* normalization                                                */
   /* ------------------------------------------------------------ */

   ML_DVector_GetDataPtr(&(ml->Amat_Normalization[level]), &scales) ;

/* watch out for this bomb */
scales = NULL;
   if ( scales != NULL ) {
      for ( i = 0; i < leng; i++ ) din_temp[i] = din[i] / scales[i];
   } else {
      scales = NULL;
      for ( i = 0; i < leng; i++ ) din_temp[i] = din[i];
   }

   /* ------------------------------------------------------------ */
   /* call MG v-cycle                                              */
   /* ------------------------------------------------------------ */

   ML_Cycle_MGV(&(ml->SingleLevel[ml->ML_finest_level]), dout, din_temp, 
                ML_ZERO, ml->comm, ML_NO_RES_NORM);

   free(din_temp);
   return 0;
}

/*****************************************************************************/
/* solve using V-cycle multigrid                                             */
/*-------------------------------------------------------------------------- */

double ML_Cycle_MGV(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not)
{
   int         i, lengc, lengf;
   double      *res,  *sol2, *rhs2, res_norm = 0., *normalscales;
   double      *rhss, *dtmp;
   ML_Operator *Amat, *Rmat;
   ML_Smoother *pre,  *post;
   ML_CSolve   *csolve;
#ifdef RAP_CHECK
   double    norm1, norm2;
#endif

#ifdef ML_ANALYSIS
   short   dummy;
   int     *cols, Nrows, j, allocated_space, ncols, lwork, info;
   double  *squareA, *vals, *eig, *work;
   char    instring[100], jobz, jobz2;
   FILE    *fp;
#endif


   Amat     = curr->Amat;
   Rmat     = curr->Rmat;
   pre      = curr->pre_smoother;
   post     = curr->post_smoother;
   csolve   = curr->csolve;
   lengf    = Amat->outvec_leng;

   /* ------------------------------------------------------------ */
   /* first do the normalization                                   */
   /* ------------------------------------------------------------ */

   rhss = (double *) malloc( lengf * sizeof(double) );
   ML_DVector_GetDataPtr(curr->Amat_Normalization, &normalscales) ;
   for ( i = 0; i < lengf; i++ ) rhss[i] = rhs[i];

#ifdef ML_ANALYSIS
   if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
   {
      fp = fopen("mlmatlab.m", "w");
      Nrows = lengf;
      fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
      allocated_space = 100;
      cols = (int    *) malloc(allocated_space*sizeof(int   ));
      vals = (double *) malloc(allocated_space*sizeof(double));
      for (i = 0; i < lengf; i++) {
         while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)==0)
         {
            allocated_space = 2*allocated_space + 1;
            free(vals); free(cols);
            cols = (int    *) malloc(allocated_space*sizeof(int   ));
            vals = (double *) malloc(allocated_space*sizeof(double));
            if (vals == NULL) {
               printf("Not enough space to get matrix row. Row length of\n");
               printf("%d was not sufficient\n",(allocated_space-1)/2);
               exit(1);
            }
         }
         for (j = 0; j < ncols; j++)
            fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
      }
      fprintf(fp, "[eigv,eig]=eig(full(A));\n");
      fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
      for (j = 0; j < Nrows; j++)
         fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,rhss[j]);
      fprintf(fp, "res=eigv'*rhs;\n");
      fprintf(fp, "plot(res)\n");
      fclose(fp);
      printf("** BEFORE pre-smoothing -- \n");
      printf("** Now you can use matlab to call the file mlmatlab.m\n");
      printf("   and it should display the residual for each eigenvalue.\n");
      printf("Press y when you are done.");
      scanf("%s", instring);
   }
#endif

   /* ------------------------------------------------------------ */
   /* smoothing or coarse solve                                    */
   /* ------------------------------------------------------------ */
   if (Rmat->to == NULL) {    /* coarsest grid */
      if ( ML_CSolve_Check( csolve ) == 1 ) {
         ML_CSolve_Apply(csolve, lengf, sol, lengf, rhss);
      } else {
         ML_Smoother_Apply(pre, lengf, sol, lengf, rhss, approx_all_zeros);
         ML_Smoother_Apply(post, lengf, sol, lengf, rhss, ML_NONZERO);
      }
      if (res_norm_or_not == ML_COMPUTE_RES_NORM) {
         res = (double *) malloc(lengf*sizeof(double));
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhss[i] - res[i];
         res_norm = sqrt(ML_gdot(lengf, res, res, comm));
         free(res);
      }
   }
   else {
      res = (double *) malloc(lengf*sizeof(double));

      /* --------------------------------------------------------- */
      /* pre-smoothing and compute residual                        */
      /* --------------------------------------------------------- */
      ML_Smoother_Apply(pre, lengf, sol, lengf, rhss, approx_all_zeros);

      if ( ( approx_all_zeros != ML_ZERO ) || 
           ( pre->smoother->ML_id != ML_EMPTY ) ) 
      {
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhss[i] - res[i];
      }
      else for ( i = 0; i < lengf; i++ ) res[i] = rhss[i];

#ifdef ML_ANALYSIS
      if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
      {
         fp = fopen("mlmatlab.m", "w");
         Nrows = lengf;
         fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
         allocated_space = 100;
         cols = (int    *) malloc(allocated_space*sizeof(int   ));
         vals = (double *) malloc(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               free(vals); free(cols);
               cols = (int    *) malloc(allocated_space*sizeof(int   ));
               vals = (double *) malloc(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
         }
         fprintf(fp, "[eigv,eig]=eig(full(A));\n");
         fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
         for (j = 0; j < Nrows; j++)
            fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,res[j]);
         fprintf(fp, "res=eigv'*rhs;\n");
         fprintf(fp, "plot(res)\n");
         fclose(fp);
         printf("** AFTER  pre-smoothing -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
     }
/*
         printf("Constructing eigensystem...\n");
         Nrows = lengf;
         squareA = malloc( Nrows * Nrows * sizeof(double) );
         for ( i = 0; i < Nrows*Nrows; i++ ) squareA[i] = 0.0;
         allocated_space = 100;
         cols = (int    *) malloc(allocated_space*sizeof(int   ));
         vals = (double *) malloc(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               free(vals); free(cols);
               cols = (int    *) malloc(allocated_space*sizeof(int   ));
               vals = (double *) malloc(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               squareA[cols[j]*Nrows+i] = vals[j];
         }
         free(cols); free(vals);
         eig = (double *) malloc( Nrows * sizeof(double) );
         work = (double *) malloc( 4 * Nrows * sizeof(double) );
         lwork = 4 * Nrows;
         jobz = 'V'; jobz2 = 'U';
         dsyev_(&jobz,&jobz2,&Nrows,squareA,&Nrows,eig,work,&lwork,&info,
                dummy,dummy);
         printf("returning from dsyev ...\n");
         if ( info > 0 )
         {
            printf("No convergence in computing the eigenvalues.\n");
         } else if ( info < 0 )
         {
            printf("Eigenvalue computation: %d-th argument has error.\n",-info);
         } else {
            fp = fopen("mlmatlab.eig", "w");
            fprintf(fp, "%d\n", Nrows);
            for (i=0; i<Nrows; i++) fprintf(fp, "%25.16e\n",eig[i]);
            for (i=0; i<Nrows*Nrows; i++)
               fprintf(fp, "%25.16e\n",squareA[i]);
            fclose(fp);
         }
         free(squareA);
         free(eig);
         free(work);
         fp = fopen("mlmatlab.m", "w");
         fprintf(fp, "res = [\n");
         for (i=0; i<Nrows; i++)
         {
            work[i] = 0.0;
            for (j=0; j<Nrows; j++)
               work[i] += ( squareA[i*Nrows+j] * rhss[j]);
            fprintf(fp, "    %25.16e \n", work[i]);
         }
         fclose(fp);
         printf("** BEFORE pre-smoothing -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
*/
#endif

      if (res_norm_or_not == ML_COMPUTE_RES_NORM)
         res_norm = sqrt(ML_gdot(lengf, res, res, comm));

      lengc = Rmat->outvec_leng;

      rhs2 = (double *) malloc(lengc*sizeof(double));
      sol2 = (double *) malloc(lengc*sizeof(double));
      for ( i = 0; i < lengc; i++ ) sol2[i] = 0.0;

      /* --------------------------------------------------------- */
      /* normalization                                             */
      /* --------------------------------------------------------- */
      ML_DVector_GetDataPtr(curr->Amat_Normalization, &normalscales) ;
      if ( normalscales != NULL )
         for ( i = 0; i < lengf; i++ ) res[i] /= normalscales[i];

      /* ------------------------------------------------------------ */
      /* transform the data from equation to grid space, do grid      */
      /* transfer and then transfer back to equation space            */
      /* ------------------------------------------------------------ */
      if ( ML_Mapper_Check(curr->eqn2grid) == 1 )
      {
         dtmp = (double *) malloc( lengf * sizeof( double ) );
         ML_Mapper_Apply(curr->eqn2grid, res, dtmp );
         for ( i = 0; i < lengf; i++ ) res[i] = dtmp[i];
         free( dtmp );
      }
      ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, res, lengc, rhs2);
      if ( ML_Mapper_Check(Rmat->to->grid2eqn) == 1 )
      {
         dtmp = (double *) malloc( lengc * sizeof( double ) );
         ML_Mapper_Apply(Rmat->to->grid2eqn, rhs2, dtmp );
         for ( i = 0; i < lengc; i++ ) rhs2[i] = dtmp[i];
         free( dtmp );
      }
      ML_DVector_GetDataPtr(Rmat->to->Amat_Normalization,&normalscales);
      if ( normalscales != NULL )
         for ( i = 0; i < lengc; i++ ) rhs2[i] = rhs2[i] * normalscales[i];

      /* --------------------------------------------------------- */
      /* process the next level and transfer back to this level    */
      /* --------------------------------------------------------- */
      ML_Cycle_MGV( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM);

      /* ------------------------------------------------------------ */
      /* transform the data from equation to grid space, do grid      */
      /* transfer and then transfer back to equation space            */
      /* ------------------------------------------------------------ */
      if ( ML_Mapper_Check(Rmat->to->eqn2grid) == 1 )
      {
         dtmp = (double *) malloc( lengc * sizeof( double ) );
         ML_Mapper_Apply(Rmat->to->eqn2grid, sol2, dtmp);
         for ( i = 0; i < lengc; i++ ) sol2[i] = dtmp[i];
         free( dtmp );
      }
      ML_Operator_ApplyAndResetBdryPts(Rmat->to->Pmat,lengc,sol2,lengf,res);
      if ( ML_Mapper_Check(curr->grid2eqn) == 1 )
      {
         dtmp = (double *) malloc( lengf * sizeof( double ) );
         ML_Mapper_Apply(curr->grid2eqn, res, dtmp);
         for ( i = 0; i < lengf; i++ ) res[i] = dtmp[i];
         free( dtmp );
      }

      /* --------------------------------------------------------- */
      /* post-smoothing                                            */
      /* --------------------------------------------------------- */
      for ( i = 0; i < lengf; i++ ) sol[i] += res[i];
#if defined(RAP_CHECK) || defined(ANALYSIS)

   /* When using RAP, the restricted residual after the coarse grid */
   /* correction should be zero.                                    */

   ML_Operator_Apply(Amat, lengf, sol, lengf, res);
   for ( i = 0; i < lengf; i++ ) res[i] = rhs[i] - res[i];

#ifdef ML_ANALYSIS
      if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
      {
         fp = fopen("mlmatlab.m", "w");
         Nrows = lengf;
         fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
         allocated_space = 100;
         cols = (int    *) malloc(allocated_space*sizeof(int   ));
         vals = (double *) malloc(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               free(vals); free(cols);
               cols = (int    *) malloc(allocated_space*sizeof(int   ));
               vals = (double *) malloc(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
         }
         fprintf(fp, "[eigv,eig]=eig(full(A));\n");
         fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
         for (j = 0; j < ncols; j++)
            fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,res[j]);
         fprintf(fp, "res=eigv'*rhs;\n");
         fprintf(fp, "plot(res)\n");
         fclose(fp);
         printf("** AFTER  coarse grid correction -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
     }
#endif

   ML_DVector_GetDataPtr(Rmat->from->Amat_Normalization,&normalscales);

   if ( normalscales != NULL )
      for ( i = 0; i < lengf; i++ ) res[i] = res[i]/normalscales[i];

   ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, res, lengc, rhs2);

   ML_DVector_GetDataPtr(Rmat->to->Amat_Normalization,&normalscales);
   if ( normalscales != NULL )
      for ( i = 0; i < lengc; i++ ) rhs2[i] = rhs2[i] * normalscales[i];

   norm1 = sqrt(ML_gdot(lengc, rhs2, rhs2, comm));
   norm2 = sqrt(ML_gdot(lengf, res, res, comm));
   if (comm->ML_mypid == 0) printf("|R r| = %e, |r| =  %e\n",norm1, norm2);

#endif
      free(sol2);
      free(rhs2);

      ML_Smoother_Apply(post, lengf, sol, lengf, rhss, ML_NONZERO);
#ifdef ML_ANALYSIS
      if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
      {
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhss[i] - res[i];
         fp = fopen("mlmatlab.m", "w");
         Nrows = lengf;
         fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
         allocated_space = 100;
         cols = (int    *) malloc(allocated_space*sizeof(int   ));
         vals = (double *) malloc(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               free(vals); free(cols);
               cols = (int    *) malloc(allocated_space*sizeof(int   ));
               vals = (double *) malloc(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
         }
         fprintf(fp, "[eigv,eig]=eig(full(A));\n");
         fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
         for (j = 0; j < ncols; j++)
            fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,res[j]);
         fprintf(fp, "res=eigv'*rhs;\n");
         fprintf(fp, "plot(res)\n");
         fclose(fp);
         printf("** AFTER  postsmoothing -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
     }
#endif
      free(res);
   }

   free(rhss);
   return(res_norm);
}

/*****************************************************************************/
/* solve using V-cycle algebraic multigrid                                   */
/*-------------------------------------------------------------------------- */

int ML_Solve_AMGV( ML *ml , double *din, double *dout)
{
   int    i, leng, level;

   /* ------------------------------------------------------------ */
   /* initially set the solution to be all 0                       */
   /* ------------------------------------------------------------ */

   level = ml->ML_finest_level;
   leng = ml->Amat[level].outvec_leng;
   for ( i = 0; i < leng; i++ ) dout[i] = 0.0;

   /* ------------------------------------------------------------ */
   /* call MG v-cycle                                              */
   /* ------------------------------------------------------------ */

   ML_Cycle_AMGV(&(ml->SingleLevel[ml->ML_finest_level]), dout, din, 
                ML_ZERO, ml->comm);

   return 0;
}

/*****************************************************************************/
/* function to perform one AMG V-cycle                                       */
/* ------------------------------------------------------------------------- */

double ML_Cycle_AMGV(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm)
{
   int         i, lengc, lengf;
   double      *res,  *sol2, *rhs2, res_norm = 0.;
   ML_Operator *Amat, *Rmat;
   ML_Smoother *pre,  *post;
   ML_CSolve   *csolve;
   static int fine_size = 0;

#ifdef RAP_CHECK
   double    norm1, norm2;
#endif

#ifdef ML_ANALYSIS
   short   dummy;
   int     *cols, Nrows, j, allocated_space, ncols, lwork, info;
   double  *squareA, *vals, *eig, *work; 
   char    instring[100], jobz, jobz2;
   FILE    *fp;
#endif

   Amat     = curr->Amat;
   Rmat     = curr->Rmat;
   pre      = curr->pre_smoother;
   post     = curr->post_smoother;
   csolve   = curr->csolve;
   lengf    = Amat->outvec_leng;

   if (fine_size == 0) fine_size = lengf;

#ifdef ML_ANALYSIS
   if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
   {
      fp = fopen("mlmatlab.m", "w");
      Nrows = lengf;
      fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
      allocated_space = 100;
      cols = (int    *) malloc(allocated_space*sizeof(int   ));
      vals = (double *) malloc(allocated_space*sizeof(double));
      for (i = 0; i < lengf; i++) {
         while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)==0) 
         {
            allocated_space = 2*allocated_space + 1;
            free(vals); free(cols);
            cols = (int    *) malloc(allocated_space*sizeof(int   ));
            vals = (double *) malloc(allocated_space*sizeof(double));
            if (vals == NULL) {
               printf("Not enough space to get matrix row. Row length of\n");
               printf("%d was not sufficient\n",(allocated_space-1)/2);
               exit(1);
            }
         }
         for (j = 0; j < ncols; j++)
            fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
      }
      free(cols);
      free(vals);
      fprintf(fp, "[eigv,eig]=eig(full(A));\n");
      fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
      for (j = 0; j < Nrows; j++)
         fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,rhs[j]);
      fprintf(fp, "res=eigv'*rhs;\n");
      fprintf(fp, "plot(res)\n");
      fclose(fp);
      printf("** BEFORE pre-smoothing -- \n");
      printf("** Now you can use matlab to call the file mlmatlab.m\n");
      printf("   and it should display the residual for each eigenvalue.\n");
      printf("Press y when you are done.");
      scanf("%s", instring);
   }
#endif

   /* ------------------------------------------------------------ */
   /* smoothing or coarse solve                                    */
   /* ------------------------------------------------------------ */

   if (Rmat->to == NULL)     /* coarsest grid */
   {
      if ( ML_CSolve_Check( csolve ) == 1 ) {
         ML_CSolve_Apply(csolve, lengf, sol, lengf, rhs);
      } else {
         ML_Smoother_Apply(pre, lengf, sol, lengf, rhs, approx_all_zeros);
         ML_Smoother_Apply(post, lengf, sol, lengf, rhs, ML_NONZERO);
      }
      if ( (lengf == fine_size) && (curr->Pmat->to == NULL)) {
         res = (double *) malloc(lengf*sizeof(double));
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhs[i] - res[i];
         res_norm = sqrt(ML_gdot(lengf, res, res, comm));
         free(res);
      }
   }
   else 
   {
      res = (double *) malloc(lengf*sizeof(double));

      /* --------------------------------------------------------- */
      /* pre-smoothing and compute residual                        */
      /* --------------------------------------------------------- */

      ML_Smoother_Apply(pre, lengf, sol, lengf, rhs, approx_all_zeros);

      if ( ( approx_all_zeros != ML_ZERO ) || 
           ( pre->smoother->ML_id != ML_EMPTY ) ) 
      {
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhs[i] - res[i];
      }
      else for ( i = 0; i < lengf; i++ ) res[i] = rhs[i];

#ifdef ML_ANALYSIS
      if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
      {
         fp = fopen("mlmatlab.m", "w");
         Nrows = lengf;
         fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
         allocated_space = 100;
         cols = (int    *) malloc(allocated_space*sizeof(int   ));
         vals = (double *) malloc(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0) 
            {
               allocated_space = 2*allocated_space + 1;
               free(vals); free(cols);
               cols = (int    *) malloc(allocated_space*sizeof(int   ));
               vals = (double *) malloc(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
         }
         fprintf(fp, "[eigv,eig]=eig(full(A));\n");
         fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
         for (j = 0; j < Nrows; j++)
            fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,res[j]);
         fprintf(fp, "res=eigv'*rhs;\n");
         fprintf(fp, "plot(res)\n");
         fclose(fp);
         printf("** AFTER  pre-smoothing -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
     }
/*
         printf("Constructing eigensystem...\n");
         Nrows = lengf;
         squareA = malloc( Nrows * Nrows * sizeof(double) );
         for ( i = 0; i < Nrows*Nrows; i++ ) squareA[i] = 0.0;
         allocated_space = 100;
         cols = (int    *) malloc(allocated_space*sizeof(int   ));
         vals = (double *) malloc(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) 
         {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols) 
                  == 0) 
            {
               allocated_space = 2*allocated_space + 1;
               free(vals); free(cols);
               cols = (int    *) malloc(allocated_space*sizeof(int   ));
               vals = (double *) malloc(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               squareA[cols[j]*Nrows+i] = vals[j];
         }
         free(cols); free(vals);
         eig = (double *) malloc( Nrows * sizeof(double) );
         work = (double *) malloc( 4 * Nrows * sizeof(double) );
         lwork = 4 * Nrows;
         jobz = 'V'; jobz2 = 'U';
         dsyev_(&jobz,&jobz2,&Nrows,squareA,&Nrows,eig,work,&lwork,&info,
                dummy,dummy);
         printf("returning from dsyev ...\n");
         if ( info > 0 ) 
         {
            printf("No convergence in computing the eigenvalues.\n");
         } 
         else if ( info < 0 )
         {
            printf("Eigenvalue computation:%d-th argument has error.\n",-info);
         } 
         else 
         {
            fp = fopen("mlmatlab.eig", "w");
            fprintf(fp, "%d\n", Nrows);
            for (i=0; i<Nrows; i++) fprintf(fp, "%25.16e\n",eig[i]);
            for (i=0; i<Nrows*Nrows; i++) 
               fprintf(fp, "%25.16e\n",squareA[i]);
            fclose(fp);
         }
         free(squareA);
         free(eig);
         free(work);
         fp = fopen("mlmatlab.m", "w");
         fprintf(fp, "res = [\n");
         for (i=0; i<Nrows; i++) 
         {
            work[i] = 0.0;
            for (j=0; j<Nrows; j++) 
               work[i] += ( squareA[i*Nrows+j] * rhs[j]); 
            fprintf(fp, "    %25.16e \n", work[i]);
         }
         fclose(fp);
         printf("** BEFORE pre-smoothing -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
*/
#endif

      if ( (lengf == fine_size) && (curr->Pmat->to == NULL))
         res_norm = sqrt(ML_gdot(lengf, res, res, comm));

      lengc = Rmat->outvec_leng;

      if ( lengc > 0 )
      {
         rhs2 = (double *) malloc(lengc*sizeof(double));
         sol2 = (double *) malloc(lengc*sizeof(double));
      }
      for ( i = 0; i < lengc; i++ ) sol2[i] = 0.0;

      /* ------------------------------------------------------------ */
      /* perform grid transfer                                        */
      /* ------------------------------------------------------------ */

      ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, res, lengc, rhs2);

      /* --------------------------------------------------------- */
      /* process the next level and transfer back to this level    */
      /* --------------------------------------------------------- */

      ML_Cycle_AMGV( Rmat->to, sol2, rhs2, ML_ZERO, comm);

      /* ------------------------------------------------------------ */
      /* transform the data from equation to grid space, do grid      */
      /* transfer and then transfer back to equation space            */
      /* ------------------------------------------------------------ */

      ML_Operator_ApplyAndResetBdryPts(Rmat->to->Pmat,lengc,sol2,lengf,res);

      /* --------------------------------------------------------- */
      /* post-smoothing                                            */
      /* --------------------------------------------------------- */
      for ( i = 0; i < lengf; i++ ) sol[i] += res[i];

#if defined(RAP2_CHECK) || defined(ANALYSIS)

      /* When using RAP, the restricted residual after the coarse grid */
      /* correction should be zero.                                    */

      printf("RAPCHECK\n");
      ML_Operator_Apply(Amat, lengf, sol, lengf, res);
      for ( i = 0; i < lengf; i++ ) res[i] = rhs[i] - res[i];

#ifdef ML_ANALYSIS
      if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
      {
         fp = fopen("mlmatlab.m", "w");
         Nrows = lengf;
         fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
         allocated_space = 100;
         cols = (int    *) malloc(allocated_space*sizeof(int   ));
         vals = (double *) malloc(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols) 
                  == 0) 
            {
               allocated_space = 2*allocated_space + 1;
               free(vals); free(cols);
               cols = (int    *) malloc(allocated_space*sizeof(int   ));
               vals = (double *) malloc(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
         }
         fprintf(fp, "[eigv,eig]=eig(full(A));\n");
         fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
         for (j = 0; j < ncols; j++)
            fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,res[j]);
         fprintf(fp, "res=eigv'*rhs;\n");
         fprintf(fp, "plot(res)\n");
         fclose(fp);
         printf("** AFTER  coarse grid correction -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
     }
#endif

      ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, res, lengc, rhs2);

      ML_DVector_GetDataPtr(Rmat->to->Amat_Normalization,&normalscales);
      if ( normalscales != NULL )
         for ( i = 0; i < lengc; i++ ) rhs2[i] = rhs2[i] * normalscales[i];

      norm1 = sqrt(ML_gdot(lengc, rhs2, rhs2, comm));
      norm2 = sqrt(ML_gdot(lengf, res, res, comm));
      if (comm->ML_mypid == 0) printf("|R r| = %e, |r| =  %e\n",norm1, norm2);

#endif
      if ( lengc > 0 ) free(sol2);
      if ( lengc > 0 ) free(rhs2);

      ML_Smoother_Apply(post, lengf, sol, lengf, rhs, ML_NONZERO);

#ifdef ML_ANALYSIS
      if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
      {
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhs[i] - res[i];
         fp = fopen("mlmatlab.m", "w");
         Nrows = lengf;
         fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
         allocated_space = 100;
         cols = (int    *) malloc(allocated_space*sizeof(int   ));
         vals = (double *) malloc(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols) 
                  == 0) 
            {
               allocated_space = 2*allocated_space + 1;
               free(vals); free(cols);
               cols = (int    *) malloc(allocated_space*sizeof(int   ));
               vals = (double *) malloc(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
         }
         free(cols); free(vals);
         fprintf(fp, "[eigv,eig]=eig(full(A));\n");
         fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
         for (j = 0; j < ncols; j++)
            fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,res[j]);
         fprintf(fp, "res=eigv'*rhs;\n");
         fprintf(fp, "plot(res)\n");
         fclose(fp);
         printf("** AFTER  postsmoothing -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
     }
#endif
      free(res);
   }
   return(res_norm);
}

/*****************************************************************************/
/*****************************************************************************/
/* SuperLU stuff                                                             */
/*****************************************************************************/
/*****************************************************************************/

#ifdef SUPERLU
#include "dsp_defs.h"
#include "util.h"
#elif DSUPERLU
#include "superlu_ddefs.h"
#endif


/*****************************************************************************/
/* Form global matrix in CSR                                                 */
/*****************************************************************************/

int ML_Gen_Amatrix_Global(ML_Matrix_DCSR *inmat, ML_Matrix_DCSR *outmat, 
                          ML_Comm *comm, int *offset) 
{
   int        nprocs, mypid, *mat_ia, *mat_ja, N_internal, N_total;
   int        i, j, k, nnz, *mat2_ia, *mat2_ja, N_external;
   int        index, *itmp, *proc_array;
   int        cur_nrows, new_nrows, cur_nnz, new_nnz;
   double     *mat_a, *mat2_a, *dtmp;

   /* ---------------------------------------------------------------- */
   /* fetch parallel machine parameters                                */
   /* ---------------------------------------------------------------- */

   nprocs = comm->ML_nprocs;
   mypid  = comm->ML_mypid;

   /* ---------------------------------------------------------------- */
   /* fetch matrix components and if it is in VBR format, first have   */
   /* to convert it to MSR format                                      */
   /* ---------------------------------------------------------------- */

   mat_ia     = inmat->mat_ia;
   mat_ja     = inmat->mat_ja;
   mat_a      = inmat->mat_a;
   N_internal = inmat->mat_n;
   nnz        = mat_ia[N_internal];
   if ( inmat->comminfo->neighbors != NULL ) { 
      N_external  = inmat->comminfo->total_rcv_length;
      N_total     = N_internal + N_external;
   } else {
      N_external  = 0;
      N_total     = N_internal;
   }

   /* ---------------------------------------------------------------- */
   /* collect the correct indices for the reconstruction of the        */
   /* global matrix                                                    */
   /* ---------------------------------------------------------------- */

   ML_memory_alloc( (void**) &proc_array, nprocs * sizeof(int), "KLA" );
   ML_memory_alloc( (void**) &itmp, nprocs * sizeof(int), "KLB" );
   for ( i = 0; i < nprocs; i++ ) proc_array[i] = 0;
   proc_array[mypid] = N_internal;
   ML_Comm_GsumVecInt(comm, proc_array, itmp, nprocs);
   for ( i = nprocs-1; i >= 1; i-- ) proc_array[i] = proc_array[i-1];
   proc_array[0] = 0;
   for ( i = 1; i < nprocs; i++ ) proc_array[i] += proc_array[i-1];
   ML_memory_free( (void**) &itmp );
   (*offset) = proc_array[mypid];
   ML_memory_alloc( (void**) &dtmp, N_total * sizeof(double), "KLC" );
   for ( i = 0; i < N_internal; i++ ) 
      dtmp[i] = (double) (proc_array[mypid] + i);
   ML_exchange_bdry(dtmp, inmat->comminfo, N_internal, comm, ML_OVERWRITE);
   ML_memory_alloc( (void**) &itmp, N_total * sizeof(int), "KLE" );
   for ( i = 0; i < N_total; i++ ) itmp[i] = (int) dtmp[i];
   ML_memory_free( (void **) &dtmp );
   ML_memory_free( (void **) &proc_array );

   /* ---------------------------------------------------------------- */
   /* communicate the sub-parts of the global matrix                   */
   /* ---------------------------------------------------------------- */

   cur_nrows = N_internal;
   new_nrows = ML_Comm_GsumInt(comm, cur_nrows);
   cur_nnz   = nnz;
   new_nnz   = ML_Comm_GsumInt(comm, cur_nnz);
   ML_memory_alloc( (void**) &mat2_ia, (new_nrows+1)*sizeof(int), "KLF");
   ML_memory_alloc( (void**) &mat2_ja, new_nnz * sizeof(int), "KLG");
   ML_memory_alloc( (void**) &mat2_a, new_nnz * sizeof(double), "KLH");
   index = 0;
   for ( i = 0; i < N_internal; i++ )
   {
      for ( j = mat_ia[i]; j < mat_ia[i+1]; j++ )
      {
         k = mat_ja[j];
         mat2_ja[index] = itmp[k];
         mat2_a[index++] = mat_a[j];  
      }
      mat2_ia[i] = mat_ia[i+1] - mat_ia[i];  
   }    
   ML_memory_free( (void **) &itmp );

   k = cur_nrows;
   ML_Comm_GappendInt(comm, mat2_ia, &k, new_nrows);
   k = index;
   ML_Comm_GappendInt(comm, mat2_ja, &k, new_nnz);
   k = index;
   ML_Comm_GappendDouble(comm, mat2_a, &k, new_nnz);

   /* ---------------------------------------------------------------- */
   /* store the incoming global matrix                                 */
   /* ---------------------------------------------------------------- */

   for ( i = 1; i < new_nrows; i++ ) mat2_ia[i] += mat2_ia[i-1];
   for ( i = new_nrows; i > 0; i-- ) mat2_ia[i] = mat2_ia[i-1];
   mat2_ia[0] = 0;

#ifdef ML_PRINT_COARSE_MAT
if ( comm->ML_mypid == 0 && new_nrows == comm->ML_nprocs)
for ( i = 0; i < new_nrows; i++ ) 
{
   for ( j = mat2_ia[i]; j < mat2_ia[i+1]; j++ ) 
      printf("A(%4d,%4d) = %e;\n", i+1, mat2_ja[j]+1, mat2_a[j]);
}
#endif

   outmat->mat_n     = new_nrows;
   outmat->mat_ia    = mat2_ia;
   outmat->mat_ja    = mat2_ja;
   outmat->mat_a     = mat2_a;

   return 0;
}

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
#elif DSUPERLU
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
int nblocks = 1, *block_list, old_upper, count, newptr, me, nnzs;
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
   row_ptr = (int *) malloc(sizeof(int)*(osize+1));
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
      cols    = (int    *) malloc(sizeof(int)*space);
      vals    = (double *) malloc(sizeof(double)*space);

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
   csr_mat = (ML_Matrix_DCSR *) malloc(sizeof(ML_Matrix_DCSR));
   csr_mat->mat_n  = osize;
   csr_mat->mat_ja = cols;
   csr_mat->mat_a  = vals;
   csr_mat->mat_ia = row_ptr;
   csr_mat->comminfo = op->getrow->pre_comm;

   /* ----------------------------------------------------------------- */
   /* form a global matrix                                              */
   /* ----------------------------------------------------------------- */

   csr2_mat = (ML_Matrix_DCSR *) malloc(sizeof(ML_Matrix_DCSR));
   ML_Gen_Amatrix_Global( csr_mat, csr2_mat, ml_handle->comm, &offset);
   free(row_ptr);
   free(cols);
   free(vals);
   free(csr_mat);

   /* Throw away some information to make it cheaper for LU. We do this   */ 
   /* by using metis to generate some blocks and factor the block matrix. */
   if (nblocks > 1) { 
      mat_ia  = csr2_mat->mat_ia;
      mat_ja  = csr2_mat->mat_ja;
      mat_val = csr2_mat->mat_a;
      nrows   = csr2_mat->mat_n;
      temp_ptr =(struct ML_CSR_MSRdata *) malloc(sizeof(struct ML_CSR_MSRdata));
      temp_ptr->rowptr = mat_ia;
      temp_ptr->columns= mat_ja;
      temp_ptr->values = mat_val;
      ML_Create(&subml, 1);
      ML_Init_Amatrix(subml, 0, nrows, nrows, (void *) temp_ptr);
      ML_Set_Amatrix_Matvec(subml, 0, CSR_matvec);
      ML_CommInfoOP_Set_neighbors(&(subml->Amat[0].getrow->pre_comm), 0,
                               NULL, ML_OVERWRITE, NULL, 0);
      ML_Operator_Set_Getrow(&(subml->Amat[0]), ML_EXTERNAL, 
                             subml->Amat[0].outvec_leng, CSR_getrows);
      ML_Gen_Blocks_Metis(subml, 0, &nblocks, &block_list);
      ML_Destroy(&subml);
      free(temp_ptr);
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
   if ( sl->csolve->func->internal == SuperLU_Solve ) reuse = 1;
   else
   {
      reuse = 0;
      sl->csolve->func->internal = SuperLU_Solve;
      sl->csolve->func->ML_id = ML_INTERNAL;
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
         SuperLU_Solve((void*)solver, 0, NULL, 0, NULL);
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
   free(csr2_mat);
#ifdef ML_TIMING
   sl->csolve->build_time = GetClock() - t0;
   ml_handle->timing->total_build_time += sl->csolve->build_time;
#endif

#elif DSUPERLU
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
int nblocks = 1, *block_list, old_upper, count, newptr, me, nnzs;
   ML *subml;

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
   row_ptr = (int *) malloc(sizeof(int)*(osize+1));
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
      cols    = (int    *) malloc(sizeof(int)*space);
      vals    = (double *) malloc(sizeof(double)*space);

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
   csr_mat = (ML_Matrix_DCSR *) malloc(sizeof(ML_Matrix_DCSR));
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

   csr2_mat = (ML_Matrix_DCSR *) malloc(sizeof(ML_Matrix_DCSR));
   ML_Gen_Amatrix_Global( csr_mat, csr2_mat, ml_handle->comm, &offset);
   free(cols);
   free(vals);
   free(row_ptr);
   free(csr_mat);

   /* Throw away some information to make it cheaper for LU. We do this   */ 
   /* by using metis to generate some blocks and factor the block matrix. */
   if (nblocks > 1) { 
      mat_ia  = csr2_mat->mat_ia;
      mat_ja  = csr2_mat->mat_ja;
      mat_val = csr2_mat->mat_a;
      nrows   = csr2_mat->mat_n;
      temp_ptr =(struct ML_CSR_MSRdata *) malloc(sizeof(struct ML_CSR_MSRdata));
      temp_ptr->rowptr = mat_ia;
      temp_ptr->columns= mat_ja;
      temp_ptr->values = mat_val;
      ML_Create(&subml, 1);
      ML_Init_Amatrix(subml, 0, nrows, nrows, (void *) temp_ptr);
      ML_CommInfoOP_Set_neighbors(&(subml->Amat[0].getrow->pre_comm), 0,
                               NULL, ML_OVERWRITE, NULL, 0);
      ML_Operator_Set_Getrow(&(subml->Amat[0]), ML_EXTERNAL, 
                             subml->Amat[0].outvec_leng, CSR_getrows);

      ML_Set_Amatrix_Matvec(subml, 0, CSR_matvec);
      ML_Gen_Blocks_Metis(subml, 0, &nblocks, &block_list);
      ML_Destroy(&subml);
      free(temp_ptr);
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
   if ( sl->csolve->func->internal == SuperLU_Solve ) reuse = 1;
   else
   {
      reuse = 0;
      sl->csolve->func->internal = SuperLU_Solve;
      sl->csolve->func->ML_id = ML_INTERNAL;
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
   free(csr2_mat);
#else
   printf("ML : SuperLU not linked.\n");
#endif

   return 0;
}

/*****************************************************************************/
/* print the total time in ML                                                */
/* ------------------------------------------------------------------------- */

int ML_Print_Timing(ML *ml)
{
#ifdef ML_TIMING
   struct ML_Timing *timing;
   double t1, t2, t3, t4, t5, t6;

   timing = ml->timing;
   t1 = ML_gsum_double(timing->precond_apply_time, ml->comm);
   t2 = ML_gsum_double(timing->total_build_time, ml->comm);
   t3 = ML_gmax_double(timing->precond_apply_time, ml->comm);
   t4 = ML_gmax_double(timing->total_build_time, ml->comm);
   t5 = ML_gmax_double(-timing->precond_apply_time, ml->comm);
   t6 = ML_gmax_double(-timing->total_build_time, ml->comm);
   if (ml->comm->ML_mypid != 0) return(1);
   t1 = t1/((double) ml->comm->ML_nprocs);
   t2 = t2/((double) ml->comm->ML_nprocs);
   t5 = - t5;
   t6 = - t6;

   printf("\nML Timing information\n\n");
   if (t1 != 0.0) printf(" Time to apply preconditioner (average) = %e\n",t1);
   if (t3 != 0.0) printf(" Time to apply preconditioner (maximum) = %e\n",t3);
   if (t5 != 0.0) printf(" Time to apply preconditioner (minimum) = %e\n",t5);
   if (t2 != 0.0) printf(" Time to build kernels        (average) = %e\n",t2);
   if (t4 != 0.0) printf(" Time to build kernels        (maximum) = %e\n",t4);
   if (t6 != 0.0) printf(" Time to build kernels        (minimum) = %e\n",t6);
#endif
   return(0);
}

/*****************************************************************************/
/* function to generate interpolation operators using finite element basis   */
/* functions                                                                 */
/*-------------------------------------------------------------------------- */

int ML_Gen_GridXsferUsingFEBasis(ML *ml, int L1, int L2, int stride)
{
   int leng, leng2;
   ML_OperatorAGX  *xsfer_op;

#ifdef ML_TIMING
   double t0;
   t0 = GetClock();
#endif

   if (ml->SingleLevel[L1].Grid->gridfcn == NULL)
      return(pr_error("ML_Gen_GridXsferUsingFEBasis: First grid is missing.\n"));
   if (ml->SingleLevel[L2].Grid->gridfcn == NULL)
      return(pr_error("ML_Gen_GridXsferUsingFEBasis: Second grid is missing.\n"));
   ML_setup_grid_xsfer_op((void*) ml->SingleLevel[L1].Grid->Grid,
                          ml->SingleLevel[L1].Grid->gridfcn,
                          (void*) ml->SingleLevel[L2].Grid->Grid,
                          ml->SingleLevel[L2].Grid->gridfcn,
                          (void **) &xsfer_op, ml->comm);
   leng = ml->SingleLevel[L1].Grid->gridfcn->USR_grid_get_nvertices(
                                      ml->SingleLevel[L1].Grid->Grid);
   ML_Operator_Set_1Levels(&(ml->Rmat[L1]), &(ml->SingleLevel[L1]), 
                           &(ml->SingleLevel[L2]));
   leng2 = xsfer_op->Nlocal_rows * stride;
   ML_Operator_Set_ApplyFuncData(&(ml->Rmat[L1]),leng*stride,leng2,
                          ML_INTERNAL, (void *) xsfer_op, 
                          xsfer_op->Nlocal_rows,
                          ML_OperatorAGX_Restrict, 1);

   ML_Operator_Set_Getrow(&(ml->Rmat[L1]), ML_INTERNAL, 
	    (xsfer_op->Nlocal_rows + xsfer_op->Nremote_rows) *stride,
			  ML_OperatorAGX_Getrows);
   ml->Rmat[L1].data_destroy = ML_Operator2AGX_Destroy;

   ML_Operator_Set_1Levels(&(ml->Pmat[L2]), &(ml->SingleLevel[L2]), 
                           &(ml->SingleLevel[L1]));
   ML_Operator_Set_ApplyFuncData(&(ml->Pmat[L2]), leng2, leng*stride, 
                      ML_INTERNAL, (void *) xsfer_op, leng, 
                      ML_OperatorAGX_Prolongate, 0);
   ML_Operator_Set_Getrow(&(ml->Pmat[L2]), ML_INTERNAL, 
            ml->Pmat[L2].outvec_leng, ML_OperatorAGX_Getcols);
   xsfer_op->AGX_stride = stride;

   ML_OperatorAGX_Gen_ComminfoOp(xsfer_op, &(ml->Rmat[L1]),
        &(ml->Pmat[L2]));

#ifdef ML_TIMING
   t0 = GetClock() - t0;
   ml->timing->total_build_time   += t0;
   t0 = t0/2;
   ml->Pmat[L2].build_time = t0;
   ml->Rmat[L1].build_time = t0;
#endif
   
   return 0;
}
int ML_Gen_Blocks_Metis(ML *ml, int level, int *nblocks, int **block_list)
{
   *block_list = (int *) ML_allocate(ml->Amat[level].outvec_leng*sizeof(int));
   if (*block_list == NULL)
      pr_error("ML_Gen_Blocks_Metis: out of space\n");

   ML_Operator_BlockPartition(&(ml->Amat[level]), ml->Amat[level].outvec_leng,
                             nblocks, *block_list, NULL, NULL, 0);
   return 0;
}
