/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/*       User Interface Functions                                            */
/* ************************************************************************* */
/* ************************************************************************* */

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ml_lapack.h"
#include "ml_struct.h"
#include "ml_agg_genP.h"
#include "ml_amg_genP.h"
#include "ml_smoother.h"
#include "ml_op_utils.h"

#ifdef ML_MPI
#include "mpi.h"
#endif
#ifdef HAVE_ML_ANASAZI
extern int ML_Anasazi_Get_SpectralNorm_Anasazi(ML_Operator * Amat,
					       int MaxIters, double Tolerance,
					       int IsProblemSymmetric,
					       int UseDiagonalScaling,
					       double * LambdaMax );
#endif


/* ************************************************************************* *
 * Structure to hold user-selected ML output level.                          *
 * ************************************************************************* */

ML_PrintControl ML_PrintLevel = {0};

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

   (*ml_ptr)->ML_finest_level   = -1;
   (*ml_ptr)->ML_coarsest_level = -1;
   (*ml_ptr)->output_level    = 10;
   (*ml_ptr)->res_output_freq = 1;
   (*ml_ptr)->tolerance       = 1.e-8;
   (*ml_ptr)->max_iterations  = 1000;
   (*ml_ptr)->MinPerProc_repartition = -1;
   (*ml_ptr)->LargestMinMaxRatio_repartition = -1.;


   ML_Comm_Create( &((*ml_ptr)->comm) );
   if (global_comm == NULL)
     global_comm = (*ml_ptr)->comm;

   if ((*ml_ptr)->comm->ML_mypid == 0 && 2 < ML_Get_PrintLevel() &&
       (Nlevels <= 0)) {
     printf("ML_Create: Warning No. of requested levels = %d\n",Nlevels);
   }
   ML_memory_check("ml_create start");

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
   ML_memory_alloc((void**) &max_eigen    ,sizeof(double)*Nlevels,"MQM");
   length = sizeof(ML_DVector) * Nlevels;
   for ( i = 0; i < Nlevels; i++ ) max_eigen[i] = 0.0;
   ML_memory_alloc((void**)&Amat_Normalization, length, "MAN");

   (*ml_ptr)->ML_num_actual_levels      = -1;
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
   (*ml_ptr)->symmetrize_matrix  = ML_FALSE;
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

   if (ml != NULL)
   {
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
   }
   return 0;
}

#ifdef GREG
int ML_Destroy2(ML **ml_ptr)
{
   int i;
   ML  *ml;

   ml = (*ml_ptr);

#ifdef ML_TIMING
   if ( ml->output_level != 0 ) ML_Print_Timing(ml);
#endif

   if (ml != NULL)
   {
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
   }
   return 0;
}
#endif

/* ************************************************************************* */
/* set symmetrize option                                                     */
/* ------------------------------------------------------------------------- */

int ML_Set_Symmetrize(ML *ml, int true_or_false)
{
  if (true_or_false == ML_TRUE) ml->symmetrize_matrix = ML_TRUE;
  else if (true_or_false == ML_FALSE) ml->symmetrize_matrix = ML_FALSE;
  else {
    printf("ML_Set_Symmetrize: Unknown option %d\n",true_or_false);
    return(0);
  }
  return(1);
}

/* ************************************************************************* */
/* set debug level                                                           */
/* ------------------------------------------------------------------------- */

int ML_Set_OutputLevel(ML *ml, int output_level)
{
  ml->output_level  = output_level;
  return(1);
}

/* ************************************************************************* */
/* Get print output level.                                                   */
/* ------------------------------------------------------------------------- */

int ML_Get_PrintLevel(void)
{
  return(ML_PrintLevel.output_level);
}

/* ************************************************************************* */
/* Set print output level.                                                   */
/* ------------------------------------------------------------------------- */

int ML_Set_PrintLevel(int print_level)
{
  ML_PrintLevel.output_level = print_level;
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

int ML_Set_Comm_Send(ML *ml, int (*send)(void*,unsigned int,int,int,USR_COMM))
{
   return ( ML_Comm_Set_SendFcn (ml->comm, send ) );
}

/* ------------------------------------------------------------------------- */

int ML_Set_Comm_Recv(ML *ml, int (*recv)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*))
{
   return ( ML_Comm_Set_RecvFcn (ml->comm, recv ) );
}

/* ------------------------------------------------------------------------- */

int ML_Set_Comm_Wait(ML *ml, int (*wait)(void*,unsigned int,int*,int*,USR_COMM,USR_REQ*))
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

int ML_Set_Grid_GetDimension(ML *ml, int nl, int (*func)(void *))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetDimension(gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetNVert(ML *ml, int nl, int (*func)(void *))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetNVert( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetNElmnt(ML *ml, int nl, int (*func)(void *))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetNElmnts( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntNVert(ML *ml, int nl, int (*func)(void *, int))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntNVert( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntVertList(ML *ml, int nl, int (*func)(void *, int, int *))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntVertList( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntGlobalNum(ML *ml, int nl, ml_big_int (*func)(void *, int))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntGlobalNum( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetVertGlobalNum(ML *ml, int nl, int (*func)(void *, int))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetVertGlobalNum( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetVertCoordinate(ML *ml, int nl, int (*func)(void *, int, double *))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetVertCoordinate( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_ComputeBasisCoef(ML *ml, int nl, int (*func)(void*,int,double*,int,double*,int*))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_ComputeBasisCoef( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntVolume(ML *ml, int nl, int (*func)(void*,int,int*,double*))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntVolumes( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntMatrix(ML *ml, int nl, int (*func)(void*,int,double**))
{
   ML_GridFunc *gf;
   gf = ml->SingleLevel[nl].Grid->gridfcn;
   ML_GridFunc_Set_GetElmntMatrix( gf, func );
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Grid_GetElmntNullSpace(ML *ml, int nl, int (*func)(void*,int,double*))
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
   ML_Operator_Set_ApplyFuncData(&(ml->Amat[level]), ilen, olen, 
                             data, olen, NULL, 0);
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Amatrix_Matvec(ML *ml, int level,
                      int (*matvec)(ML_Operator *, int, double *, int, double *))
{
   ML_Operator *matrix;
   matrix = &(ml->Amat[level]);

   return(ML_Operator_Set_ApplyFunc(matrix,matvec));
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
        int (*getrow)(ML_Operator *, int , int* , int , int*, double* , int*),
	int (*comm  )(double *vec, void *data), int comm_vec_leng )
{
   ML_Operator *Amat;
   int         Nghost;

   Amat = &(ml->Amat[nl]);

   ML_Operator_Set_Getrow(Amat, Amat->outvec_leng, getrow);
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
   return 0;

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
   ML_Operator_Set_ApplyFuncData(&(ml->Rmat[level]), ilen, olen, 
                             data, olen, NULL, 0);
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Restrictor_Matvec( ML *ml , int from_level,
	int (*func) (ML_Operator *, int, double *, int, double *))
{
   ML_Operator *matrix;
   matrix = &(ml->Rmat[from_level]);

   return(ML_Operator_Set_ApplyFunc(matrix,func));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Restrictor_Getrow(ML *ml, int nl,
        int (*getrow)(ML_Operator *, int , int* , int , int*, double* , int*),
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

   return(ML_Operator_Set_Getrow(Rmat, Rmat->outvec_leng, getrow));
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
   ML_Operator_Set_ApplyFuncData(&(ml->Pmat[level]), ilen, olen, 
                             data, olen, NULL, 0);
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Set_Prolongator_Matvec( ML *ml , int to_level,
	int (*func) (ML_Operator *, int, double *, int, double *))
{
   ML_Operator *matrix;
   matrix = &(ml->Pmat[to_level]);

   return(ML_Operator_Set_ApplyFunc(matrix, func));
}

/* ------------------------------------------------------------------------- */

int ML_Set_Prolongator_Getrow(ML *ml, int nl,
        int (*getrow)(ML_Operator *, int , int* , int , int*, double* , int*),
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

   return(ML_Operator_Set_Getrow(Pmat, Pmat->outvec_leng, getrow));
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
                     int (*func)(ML_Smoother *, int, double *, int, double *),
                     char *str)
{
   char temp[80], *tptr = NULL;
   if (str != NULL) { sprintf(temp,"%s%d",str,nl); tptr = temp; }

   if (nl == ML_ALL_LEVELS) {
      printf("ML_Set_Smoother: ML_ALL_LEVELS not allowed\n");
      return 1;
   }
   if (nl < 0) {
      printf("ML_Set_Smoother: cannot set smoother on level %d\n",nl);
      return 1;
   }
   if (pre_or_post == ML_PRESMOOTHER) {
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), data,
			     func, 1,(double) ML_DEFAULT, tptr));
   }
   else if (pre_or_post == ML_POSTSMOOTHER)
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), data,
			     func, 1, (double) ML_DEFAULT, tptr));
   else if (pre_or_post == ML_BOTH)  {
      ML_Smoother_Set(&(ml->pre_smoother[nl]), data,
                      func, 1,(double) ML_DEFAULT, tptr);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), data,
                             func, 1, (double) ML_DEFAULT, tptr));
   }
   else return(pr_error("ML_Set_Smoother: unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the damped Jacobi smoother                                       */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_Jacobi( ML *ml , int nl, int pre_or_post, int ntimes,
        double omega)
{
   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int start_level, end_level, i, status = 1;
   char str[80];

   if (nl == ML_ALL_LEVELS) { start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_Smoother_Jacobi: cannot set smoother on level %d\n",start_level);
      return 1;
   }

   fun = ML_Smoother_Jacobi;
   if (omega == ML_DEFAULT) omega = .5;

   if (pre_or_post == ML_PRESMOOTHER) {
      for (i = start_level; i <= end_level; i++) {
       sprintf(str,"Jac_pre%d",i);
       status = ML_Smoother_Set(&(ml->pre_smoother[i]), NULL,
                             fun, ntimes, omega, str);
      }
   }
   else if (pre_or_post == ML_POSTSMOOTHER)
      for (i = start_level; i <= end_level; i++) {
       sprintf(str,"Jac_post%d",i);
       status = ML_Smoother_Set(&(ml->post_smoother[i]), NULL,
                             fun, ntimes,omega, str);
      }
   else if (pre_or_post == ML_BOTH) {
      for (i = start_level; i <= end_level; i++) {
        sprintf(str,"Jac_pre%d",i);
        status = ML_Smoother_Set(&(ml->pre_smoother[i]), NULL,
                             fun, ntimes, omega, str);
        sprintf(str,"Jac_post%d",i);
        status = ML_Smoother_Set(&(ml->post_smoother[i]), NULL,
                             fun, ntimes,omega, str);
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
   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int start_level, end_level, i, status = 1;
   char str[80];

   if (nl == ML_ALL_LEVELS) { start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_Smoother_GaussSeidel: cannot set smoother on level %d\n",start_level);
      return 1;
   }

   fun = ML_Smoother_GaussSeidel;

   if (pre_or_post == ML_PRESMOOTHER) {
      for (i = start_level; i <= end_level; i++) {
             sprintf(str,"GS_pre%d",i);
             status = ML_Smoother_Set(&(ml->pre_smoother[i]), NULL,
                                      fun, ntimes, omega, str);
      }
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      for (i = start_level; i <= end_level; i++) {
             sprintf(str,"GS_post%d",i);
             status = ML_Smoother_Set(&(ml->post_smoother[i]),NULL,
                             fun, ntimes, omega, str);
      }
   }
   else if (pre_or_post == ML_BOTH) {
      for (i = start_level; i <= end_level; i++) {
             sprintf(str,"GS_pre%d",i);
             status = ML_Smoother_Set(&(ml->pre_smoother[i]), NULL,
                                      fun, ntimes, omega, str);
             sprintf(str,"GS_post%d",i);
             status = ML_Smoother_Set(&(ml->post_smoother[i]), NULL,
                             fun, ntimes, omega, str);
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
   int         start_level, end_level, i, j, status = 0, Nrows, count;
   int         *bindx;
   double      *nums, **sgs_nums = NULL, *num2;
   double      *val = NULL, temp_omega;
   ML_Operator *Amat;
   int         (*fun)(ML_Smoother *, int, double *, int, double *);
   void        (*fun2)(void *) = NULL;
   struct ML_CSR_MSRdata *ptr = NULL;
   char        str[80];

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

      if (Amat->getrow->func_ptr == MSR_getrows){
         ptr   = (struct ML_CSR_MSRdata *) Amat->data;
         val   = ptr->values;
         bindx = ptr->columns;
      }

#ifdef AZTEC
      else AZ_get_MSR_arrays(Amat, &bindx, &val);
#endif
      if (val != NULL) {
         fun = ML_Smoother_MSR_SGS;
	 if (omega != 1.0) {
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
	 else {
	   fun = ML_Smoother_MSR_SGSnodamping;
	 }
      }

      if (pre_or_post == ML_PRESMOOTHER) {
         sprintf(str,"SGS_pre%d",i);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), 
                                  sgs_nums, fun, ntimes, omega, str);
	 ml->pre_smoother[i].data_destroy = fun2;
      }
      else if (pre_or_post == ML_POSTSMOOTHER) {
         sprintf(str,"SGS_post%d",i);
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), 
                                  sgs_nums, fun, ntimes, omega, str);
	 ml->post_smoother[i].data_destroy = fun2;
      }
      else if (pre_or_post == ML_BOTH) {
         sprintf(str,"SGS_pre%d",i);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), 
                                  sgs_nums, fun, ntimes, omega, str);
         sprintf(str,"SGS_post%d",i);
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), 
                                  sgs_nums, fun, ntimes, omega, str);
	 ml->post_smoother[i].data_destroy = fun2;
      }
      else return(pr_error("Print unknown pre_or_post choice\n"));
      fun2 = NULL;
      sgs_nums = NULL;
   }
   return(status);
}

/* ------------------------------------------------------------------------- */
/* generate the sequential symmetric Gauss Seidel smoother                   */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_SymGaussSeidelSequential(ML *ml , int nl, int pre_or_post,
                                             int ntimes, double omega)
{
   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int start_level, end_level, i, status;
   char str[80];

   if (nl == ML_ALL_LEVELS) { start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0)
   {
      printf("ML_Gen_Smoother_SymGaussSeidelSequential: cannot set smoother ");
      printf("on level %d\n", start_level);
      return 1;
   }

   fun = ML_Smoother_SGSSequential;
   if (omega == ML_DEFAULT) omega = 1.;

   for (i = start_level; i <= end_level; i++)
   {
      if (pre_or_post == ML_PRESMOOTHER)
      {
         sprintf(str,"SGS_pre%d",i);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]),  NULL,
                                  fun, ntimes, omega, str);
      }
      else if (pre_or_post == ML_POSTSMOOTHER)
      {
         sprintf(str,"SGS_post%d",i);
         status = ML_Smoother_Set(&(ml->post_smoother[i]),  NULL,
                                  fun, ntimes, omega, str);
      }
      else if (pre_or_post == ML_BOTH)
      {
         sprintf(str,"SGS_pre%d",i);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]),  NULL,
                                  fun, ntimes, omega, str);
         sprintf(str,"SGS_post%d",i);
         status = ML_Smoother_Set(&(ml->post_smoother[i]),  NULL,
                                  fun, ntimes, omega, str);
      }
      else return(pr_error("ML_Gen_SGSSequential: unknown pre_or_post choice\n"));
   }
   return(status);
}

/* ------------------------------------------------------------------------- */
/* generate some other Gauss Seidel smoother (Ray, what is it ?)             */
/* ------------------------------------------------------------------------- */

int ML_Gen_SmootherGSextra( ML *ml , int nl, int pre_or_post, int ntimes,
                            double omega, int Nextra, int extra[])
{
   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int start_level, end_level, i, status = 0, Nrows, count;
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

   if (Amat->getrow->func_ptr == MSR_getrows){
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


   if (pre_or_post == ML_PRESMOOTHER)
      for (i = start_level; i <= end_level; i++) {
         status = ML_Smoother_Set(&(ml->pre_smoother[i]),  sgs_nums,
				  fun, ntimes, omega, NULL);
	 ml->pre_smoother[i].data_destroy = fun2;
   }
   else if (pre_or_post == ML_POSTSMOOTHER)
      for (i = start_level; i <= end_level; i++) {
	 status = ML_Smoother_Set(&(ml->post_smoother[i]),  sgs_nums,
                             fun, ntimes, omega, NULL);
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
   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int start_level, end_level, i, status = 0, *ordering;
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

   if (pre_or_post == ML_PRESMOOTHER)
      for (i = start_level; i <= end_level; i++) {
#ifdef ML_TIMING
         t0 = GetClock();
#endif
         ML_Smoother_Gen_Ordering(&(ml->Amat[i]), &ordering);
         ml->pre_smoother[i].data_destroy = ML_Smoother_Clean_OrderedSGS;
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), 
                         (void *) ordering, fun, ntimes, omega, NULL);
#ifdef ML_TIMING
         ml->pre_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->pre_smoother[i].build_time;
#endif
      }
   else if (pre_or_post == ML_POSTSMOOTHER)
      for (i = start_level; i <= end_level; i++) {
#ifdef ML_TIMING
         t0 = GetClock();
#endif
         ML_Smoother_Gen_Ordering(&(ml->Amat[i]), &ordering);
         ml->post_smoother[i].data_destroy = ML_Smoother_Clean_OrderedSGS;
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), 
                         (void *) ordering, fun, ntimes, omega, NULL);
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
   int            (*fun)(ML_Smoother *, int, double *, int, double *);
   ML_Sm_BGS_Data *data;
   int            start_level, end_level, i, status = 1;
   char           str[80];
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

   if (pre_or_post == ML_PRESMOOTHER) {
      for (i = start_level; i <= end_level; i++) {
         ML_Smoother_Create_BGS_Data(&data);
	 ML_Smoother_Gen_BGSFacts(&data, &(ml->Amat[i]), blocksize);
	 ml->pre_smoother[i].data_destroy = ML_Smoother_Clean_BGS_Data;
         sprintf(str,"BGS_pre%d",i);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), 
		                 (void *) data, fun, ntimes, omega, str);
#ifdef ML_TIMING
         ml->pre_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->pre_smoother[i].build_time;
#endif
      }
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      for (i = start_level; i <= end_level; i++) {
         ML_Smoother_Create_BGS_Data(&data);
	 ML_Smoother_Gen_BGSFacts(&data, &(ml->Amat[i]), blocksize);
	 ml->post_smoother[i].data_destroy = ML_Smoother_Clean_BGS_Data;
         sprintf(str,"BGS_post%d",i);
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), 
			      (void *) data, fun, ntimes, omega, str);
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
         sprintf(str,"BGS_pre%d",i);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), 
		                 (void *) data, fun, ntimes, omega, str);
#ifdef ML_TIMING
         ml->pre_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->pre_smoother[i].build_time;
#endif
         sprintf(str,"BGS_post%d",i);
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), 
			      (void *) data, fun, ntimes, omega, str);
	 ml->post_smoother[i].data_destroy = ML_Smoother_Clean_BGS_Data;
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
   int            (*fun)(ML_Smoother *, int, double *, int, double *);
   double         myomega;
   ML_Sm_BGS_Data *data;
   char           str[80];

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

   if (pre_or_post == ML_PRESMOOTHER) {
      sprintf(str,"VBJ_pre%d",nl);
      ml->pre_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                        (void *) data, fun, ntimes, myomega, str));
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      sprintf(str,"VBJ_post%d",nl);
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, myomega, str));
   }
   else if (pre_or_post == ML_BOTH) {
      sprintf(str,"VBJ_pre%d",nl);
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                        (void *) data, fun, ntimes, myomega, str);
      sprintf(str,"VBJ_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, myomega, str));
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the variable block Gauss Seidel smoother                         */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockSymGaussSeidel( ML *ml , int nl, int pre_or_post,
                      int ntimes, double omega, int nBlocks, int *blockIndices)
{
   int            (*fun)(ML_Smoother *, int, double *, int, double *);
   ML_Sm_BGS_Data *data;
   char str[80];

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

   if (pre_or_post == ML_PRESMOOTHER) {
      ml->pre_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      sprintf(str,"VBSGS_pre%d",nl);
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                        (void *) data, fun, ntimes, omega, str));
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      sprintf(str,"VBSGS_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, omega, str));
   }
   else if (pre_or_post == ML_BOTH) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      sprintf(str,"VBSGS_pre%d",nl);
      ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                        (void *) data, fun, ntimes, omega, str);
      sprintf(str,"VBSGS_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, omega, str));
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the variable block Gauss Seidel smoother (sequential)            */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockSymGaussSeidelSequential( ML *ml , int nl,
     int pre_or_post, int ntimes, double omega, int nBlocks, int *blockIndices)
{
   int            (*fun)(ML_Smoother *, int, double *, int, double *);
   ML_Sm_BGS_Data *data;
   char str[80];

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

   if (pre_or_post == ML_PRESMOOTHER) {
      ml->pre_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      sprintf(str,"VBSGSS_pre%d",nl);
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                             (void *) data, fun, ntimes, omega, str));
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      sprintf(str,"VBSGSS_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, omega, str));
   }
   if (pre_or_post == ML_BOTH) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      sprintf(str,"VBSGSS_pre%d",nl);
      ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                             (void *) data, fun, ntimes, omega, str);
      sprintf(str,"VBSGSS_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, omega, str));
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the variable block Jacobi smoother with Krylov(not debugged yet) */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_VBlockKrylovJacobi( ML *ml , int nl, int pre_or_post,
                 int ntimes, double omega, int nBlocks, int *blockIndices)
{
   int            (*fun)(ML_Smoother *, int, double *, int, double *);
   double         myomega;
   ML_Sm_BGS_Data *data;
   char           str[80];

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

   if (pre_or_post == ML_PRESMOOTHER) {
      ml->pre_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      sprintf(str,"VBKJ_pre%d",nl);
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                             (void *) data, fun, ntimes, myomega,str));
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      sprintf(str,"VBKJ_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, myomega, str));
   }
   else if (pre_or_post == ML_BOTH) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_BGS_Data;
      sprintf(str,"VBKJ_pre%d",nl);
      ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                             (void *) data, fun, ntimes, myomega,str);
      sprintf(str,"VBKJ_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, myomega, str));
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

/* ------------------------------------------------------------------------- */
/* generate the overlapped domain decomposition with ILUT smoother           */
/* ------------------------------------------------------------------------- */

int ML_Gen_Smoother_OverlappedDDILUT( ML *ml , int nl, int pre_or_post )
{
   int             (*fun)(ML_Smoother *, int, double *, int, double *);
   int             total_recv_leng, *recv_lengths, *int_buf, *map, *map2;
   int             offset;
   double          *dble_buf;
   ML_Sm_ILUT_Data *data;
   ML_Operator     *Amat;
   ML_Comm         *comm;
   char            str[80];

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
   if ( map  != NULL ) ML_free(map);
   if ( map2 != NULL ) ML_free(map2);
   if ( int_buf != NULL ) ML_free(int_buf);
   if ( dble_buf != NULL ) ML_free(dble_buf);
   if ( recv_lengths != NULL ) ML_free(recv_lengths);

   /* ---------------------------------------------------------------- */
   /* set it up as smoother                                            */
   /* ---------------------------------------------------------------- */

   if (pre_or_post == ML_PRESMOOTHER) {
      ml->pre_smoother[nl].data_destroy = ML_Smoother_Destroy_ILUT_Data;
      sprintf(str,"ODDILUT_pre%d",nl);
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                        (void *) data, fun, 1, 0.0,str));
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_ILUT_Data;
      sprintf(str,"ODDILUT_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, 1, 0.0,str));
   }
   else if (pre_or_post == ML_BOTH) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_ILUT_Data;
      sprintf(str,"ODDILUT_pre%d",nl);
      ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                        (void *) data, fun, 1, 0.0,str);
      sprintf(str,"ODDILUT_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, 1, 0.0,str));
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
}

#ifdef out
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
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                        (void *) data, fun, ntimes, 0.0, str));
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_Schwarz_Data;
      sprintf(str,"VBASz_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, 0.0, str));
   }
   else if (pre_or_post == ML_BOTH) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_Schwarz_Data;
      sprintf(str,"VBASz_pre%d",nl);
      ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                      (void *) data, fun, ntimes, 0.0, str);
      sprintf(str,"VBASz_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, 0.0, str));
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
      return(ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                        (void *) data, fun, ntimes, 0.0, str));
   }
   else if (pre_or_post == ML_POSTSMOOTHER) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_Schwarz_Data;
      sprintf(str,"VBMSz_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, 0.0, str));
   }
   else if (pre_or_post == ML_BOTH) {
      ml->post_smoother[nl].data_destroy = ML_Smoother_Destroy_Schwarz_Data;
      sprintf(str,"VBMSz_pre%d",nl);
      ML_Smoother_Set(&(ml->pre_smoother[nl]), 
                        (void *) data, fun, ntimes, 0.0, str);
      sprintf(str,"VBMSz_post%d",nl);
      return(ML_Smoother_Set(&(ml->post_smoother[nl]), 
                             (void *) data, fun, ntimes, 0.0, str));
   }
   else return(pr_error("Print unknown pre_or_post choice\n"));
}
#endif


#ifdef	MB_MODIF

/* ------------------------------------------------------------------------- */
/* generate the MLS smoother */
/* ------------------------------------------------------------------------- */

#include "ml_mls.h"

int ML_MLS_Setup_Coef(void *sm, int deg, int symmetrize)
{ /*
   * Preset the coefficients for MLS smoothing on current level to
   * (rho/2)* ( 1-cos((2Pi*k)/(2*deg+1)) )
   *
   * Returns: 0 on success.
   */
   const int    nSample=20000;
   double       gridStep, rho, rho2, ddeg, aux0, aux1;
   double       aux_om, om_loc[MLS_MAX_DEG], coord, samplej;
   const double pi=4.e0 * atan(1.e0); /* 3.141592653589793115998e0; */
   int          i, j, nGrid;

   /* Get all the pointers */

   ML_Smoother *smooth_ptr = (ML_Smoother *) sm;
   ML_Operator *Amat = smooth_ptr->my_level->Amat;
   struct MLSthing *widget = (struct MLSthing *) smooth_ptr->smoother->data;

   if (deg > MLS_MAX_DEG) {
       return (pr_error("*** value of deg larger than MLS_MAX_DEG !\n"));
   }

   ML_Gimmie_Eigenvalues(Amat, ML_DIAGSCALE, ML_SYMMETRIC, symmetrize);
   rho = Amat->lambda_max;

   /* Boost the largest eigenvalue by a fudge factor  */
   /* (mlsOver is now set to 1.1 in the code[1/7/02]).*/
   /* We know that our eigenvalue estimate is a lower */
   /* bound. Undershooting the largest eigenvalue for */
   /* a polynomial method is a bit risky.             */

   rho *= widget->mlsOver;

   for (i=0; i<MLS_MAX_DEG; i++) { widget->mlsOm[i] = 0.e0; om_loc[i] = 0.e0; }

   /* compute a set of polynomial coefficients. These */
   /* almost look like Chebyshev coefficients. They   */
   /* correspond to                                   */
   /*        2/(rho * ( 1 - cos(2 pi (j+1)/(2n + 1))) */
   /* By the way, this polynomial is done via products*/
   /* I think? that it corresponds to                 */
   /*    (I - om_loc[0] A) * (I - om_loc[1] A) * ...  */

   ddeg = (double)deg;
   aux1 = 1.e0/(2.e0 * ddeg + 1.e0);

   for (j=0; j<deg; j++) {
	   aux0 = 2.e0 * (double)(j+1) * pi;
	   aux_om = rho/2.e0 * (1.e0 - cos(aux0 * aux1));
	   om_loc[j] = 1.e0/aux_om;
   }

   /* Compute the coefficients of the polynomial      */
   /*    (I - om_loc[0] A) * (I - om_loc[1] A) * ...  */

   widget->mlsCf[0] =  om_loc[0] + om_loc[1] + om_loc[2]+om_loc[3] + om_loc[4];
   widget->mlsCf[1] = -(om_loc[0]*om_loc[1]   + om_loc[0]*om_loc[2]
	              + om_loc[0]*om_loc[3]   + om_loc[0]*om_loc[4]
	              + om_loc[1]*om_loc[2]   + om_loc[1]*om_loc[3]
	              + om_loc[1]*om_loc[4]   + om_loc[2]*om_loc[3]
	              + om_loc[2]*om_loc[4]   + om_loc[3]*om_loc[4]);
   widget->mlsCf[2] = (om_loc[0]*om_loc[1]*om_loc[2]
	              + om_loc[0]*om_loc[1]*om_loc[3]
	              + om_loc[0]*om_loc[1]*om_loc[4]
	              + om_loc[0]*om_loc[2]*om_loc[3]
	              + om_loc[0]*om_loc[2]*om_loc[4]
	              + om_loc[0]*om_loc[3]*om_loc[4]
	              + om_loc[1]*om_loc[2]*om_loc[3]
	              + om_loc[1]*om_loc[2]*om_loc[4]
	              + om_loc[1]*om_loc[3]*om_loc[4]
	              + om_loc[2]*om_loc[3]*om_loc[4]);
   widget->mlsCf[3] = -(om_loc[0]*om_loc[1]*om_loc[2]*om_loc[3]
                      + om_loc[0]*om_loc[1]*om_loc[2]*om_loc[4]
                      + om_loc[0]*om_loc[1]*om_loc[3]*om_loc[4]
                      + om_loc[0]*om_loc[2]*om_loc[3]*om_loc[4]
                      + om_loc[1]*om_loc[2]*om_loc[3]*om_loc[4]);
   widget->mlsCf[4] = om_loc[0]*om_loc[1]*om_loc[2]*om_loc[3]*om_loc[4];

   /* Here I believe that we are trying to estimate the largest */
   /* to estimate the largest eigenvalue of p(A)*p(A)*A where   */
   /*   p(A) = (I - om_loc[0] A) * (I - om_loc[1] A) * ...      */
   /* Normally, deg = 1 so we just go the hardwired stuff in the*/
   /* else clause.                                              */

   if (deg > 1) {
     gridStep = rho/(double)nSample;
     nGrid    = ML_min((int)rint(rho/gridStep)+1, nSample);

     rho2  = 0.e0;
     for (j=0; j<nGrid; j++)  {
       coord   = (double)(j+1) * gridStep;
       samplej = 1.e0 - om_loc[0] * coord;
       for(i=1; i<deg; i++) {
	 samplej *= ( 1.e0 -  om_loc[i] * coord );
       }
       samplej *= samplej * coord;
       if (samplej > rho2) {rho2 = samplej; }
     }
   }
   else {
       rho2 = 4.0e0/(27.e0 * om_loc[0]);
   }

   /* Boost the eigenvalue estimate for p(A)*p(A)*A */
   /* and store coefficients for later application. */

   if (deg < 2) {
           widget->mlsBoost = 1.5e0;
	   widget->mlsBoost = 1.019e0;
   } else {
           widget->mlsBoost = 1.600e0;
	   widget->mlsBoost = 1.025e0;
   }
   rho2 *= widget->mlsBoost;
   widget->mlsOm2 = 2.0e0/rho2;
   for (i=0; i<deg; i++) widget->mlsOm[i] = om_loc[i];

   return 0;

}
int ML_Gen_Smoother_BlockDiagScaledCheby(ML *ml, int nl, int pre_or_post,
					 double eig_ratio, int deg,
					 int nBlocks, int *blockIndices)
{
  double temp;
  struct MLSthing *widget;
  ML_Operator *Amat;
  ML_Operator *blockMat;


  if (deg < 0) {
    printf("ML_Gen_Smoother_BlockDiagScaledCheby: deg < 0 not allowed\n");
    return 1;
  }
  if (nl == ML_ALL_LEVELS) {
    printf("ML_Gen_Smoother_BlockDiagScaledCheby: ML_ALL_LEVELS not supported\n");
    return 1;
  }
  Amat = &(ml->Amat[nl]);

  /* put in a bogus value for lambda so that ML_Gen_Smoother_MLS() */
  /* does not try to compute the eigenvalues of A. What we         */
  /* actually need are the eigenvalues of  Dinv A. These will be   */
  /* computed later in this function.                              */

  temp = Amat->lambda_max;
  Amat->lambda_max = 1.;
  ML_Gen_Smoother_MLS(ml, nl, pre_or_post, eig_ratio, deg);
  Amat->lambda_max = temp;

  if (pre_or_post != ML_POSTSMOOTHER) {
    widget = (struct MLSthing *) ml->pre_smoother[nl].smoother->data;
  }
  else {
    widget = (struct MLSthing *) ml->post_smoother[nl].smoother->data;
  }

  ML_Gen_BlockScaledMatrix_with_Eigenvalues(Amat, nBlocks, blockIndices,
					    &blockMat, widget);
  return 0;

}

int ML_Gen_BlockScaledMatrix_with_Eigenvalues(ML_Operator *Amat,
					      int nBlocks, 
					      int *blockIndices,
					      ML_Operator **blockMat,
					      struct MLSthing *widget) {
  ML_Sm_BGS_Data *data;

  widget->unscaled_matrix = Amat;

  ML_Smoother_Create_BGS_Data(&data);
  if (blockIndices != NULL) 
    ML_Smoother_Gen_VBGSFacts(&data, Amat, nBlocks, blockIndices);
  else
    ML_Smoother_Gen_BGSFacts(&data, Amat, Amat->num_PDEs);

  widget->block_scaling   = data;

  /* create a block matrix wrapper that will correspond to Dinv A */

  *blockMat = ML_Operator_Create(Amat->comm);
  ML_Operator_Set_ApplyFuncData(*blockMat,Amat->invec_leng, Amat->outvec_leng, 
				widget,Amat->outvec_leng, NULL,0);
  ML_Operator_Set_ApplyFunc (*blockMat, ML_BlockScaledApply);

  widget->scaled_matrix = *blockMat;

  ML_Gimmie_Eigenvalues(*blockMat, ML_NO_SCALE, ML_NONSYMM, ML_NO_SYMMETRIZE);

  return 0;

}

/* Set Amat->lambda_max and Amat->lambda_min by estimating the */
/* eigenvalues of A or D^-1 A                                  */
int ML_Gimmie_Eigenvalues(ML_Operator *Amat, int scale_by_diag,
			  int matrix_is_nonsymmetric, int symmetrize_matrix)
{
  ML_Krylov   *kdata;
  ML_Operator *t2, *t3;


  if ((Amat->lambda_max < -666.) && (Amat->lambda_max > -667)) {
     kdata = ML_Krylov_Create( Amat->comm );
     if (scale_by_diag == 0) ML_Krylov_Set_DiagScaling_Eig(kdata, 0);
     if (matrix_is_nonsymmetric && (symmetrize_matrix == 0))
#ifdef HAVE_ML_ANASAZI
       /* Ray, this is how the function is organized:
	  - 10 is the maximum number of iterations
	  - 1e10 a tolerance
	  - ML_FALSE because the problem is non-symmetric
	  - ML_TRUE because you want diagonal scaling
	  (the function is in ml_anasazi.cpp)
       */
       ML_Anasazi_Get_SpectralNorm_Anasazi(Amat,10,1e-10,
					   ML_FALSE, ML_TRUE,
					   &(Amat->lambda_max) );
#else
       ML_Krylov_Set_ComputeNonSymEigenvalues( kdata );
#endif
     else 
       ML_Krylov_Set_ComputeEigenvalues( kdata );

     ML_Krylov_Set_PrintFreq( kdata, 0 );

     if (symmetrize_matrix == ML_TRUE) {
       t2 = ML_Operator_Create(Amat->comm);
       ML_Operator_Transpose_byrow(Amat,t2);
       t3 = ML_Operator_Create(Amat->comm);
       ML_Operator_Add(Amat,t2,t3, ML_CSR_MATRIX,1.);
       ML_Krylov_Set_Amatrix(kdata, t3);
     }
     else  ML_Krylov_Set_Amatrix(kdata, Amat);
     ML_Krylov_Solve(kdata, Amat->outvec_leng, NULL, NULL);
     Amat->lambda_max = ML_Krylov_Get_MaxEigenvalue(kdata);
     Amat->lambda_min = kdata->ML_eigen_min;

     ML_Krylov_Destroy( &kdata );
     if (symmetrize_matrix == ML_TRUE) {
       if (t3 != NULL) ML_Operator_Destroy(&t3);
       if (t2 != NULL) ML_Operator_Destroy(&t2);
     }
  }
  return 0;
} 

int ML_Gen_Smoother_MLS(ML *ml, int nl, int pre_or_post,
			double eig_ratio, int deg)
{
   int              start_level, end_level, i, j, errCode=0;
   struct MLSthing *widget;
   ML_Operator     *Amat;
   double          *tdiag;
   char             str[80];
   int                (*fun)(ML_Smoother *, int, double *, int, double *);
   int iii, degree;
   int ntimes = 1;

#ifdef ML_TIMING
   double         t0;
#endif

   if (nl == ML_ALL_LEVELS) {
	  start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}

   if (start_level < 0) {
      printf("ML_Gen_Smoother_MLS: cannot set smoother on level %d\n",
		      start_level);
      return 1;
   }


   fun = ML_Smoother_MLS_Apply;
   iii = 0;
   degree = 1;


   /*
   if (ml->comm->ML_mypid == 0) {
     printf("Enter -k to do a kth degree MLS or k (where k > 0)\n");
     printf("to do kth degree Cheby over high frequencies\n");
     scanf("%d",&iii);
     iii = deg;
   }
   ML_gsum_scalar_int(&jjj, &iii, ml->comm);
   */
   iii = deg;
   if (iii >= 0) {
     fun = ML_Cheby;
     degree = iii;
   }
   else degree = -iii;



   for (i = start_level; i <= end_level; i++) {
#ifdef ML_TIMING
     t0 = GetClock();
#endif
     Amat = &(ml->Amat[i]);
     ML_Gimmie_Eigenvalues(Amat, ML_DIAGSCALE, ML_SYMMETRIC, 
			   ml->symmetrize_matrix);

     /* To avoid division by zero problem. */
     if (Amat->diagonal != NULL)
     {
        ML_Operator_Get_Diag(Amat, Amat->outvec_leng, &tdiag);
        for (j=0; j<Amat->outvec_leng; j++)
           if (tdiag[j] == 0) tdiag[j] = 1.0;
     }


     if (Amat->matvec->func_ptr != NULL) {
         widget = ML_Smoother_Create_MLS();
	 widget->mlsDeg   = degree;
	 widget->eig_ratio = eig_ratio;
	 if (pre_or_post == ML_PRESMOOTHER) {
	   ml->pre_smoother[i].data_destroy = ML_Smoother_Destroy_MLS;
           sprintf(str,"MLS_pre%d",i);
           errCode=ML_Smoother_Set(&(ml->pre_smoother[i]), 
				  (void *) widget, fun, 
				   ntimes, 0.0, str);
	   /* set up the values needed for MLS  */
	   if (fun == ML_Smoother_MLS_Apply) {
	     if (ML_MLS_Setup_Coef(&(ml->pre_smoother[i]), 1,ml->symmetrize_matrix)) {
	       return pr_error("*** MLS setup failed!\n");
	     }
	   }
	 }
	 else if (pre_or_post == ML_POSTSMOOTHER) {
	   ml->post_smoother[i].data_destroy = ML_Smoother_Destroy_MLS;
	   sprintf(str,"MLS_post%d",i);
	   errCode=ML_Smoother_Set(&(ml->post_smoother[i]), 
				  (void *) widget, fun, 
				  ntimes, 0.0, str);
	   /* set up the values needed for MLS  */
	   if (fun == ML_Smoother_MLS_Apply) {
	     if (ML_MLS_Setup_Coef(&(ml->post_smoother[i]), 1,ml->symmetrize_matrix)) {
	       return pr_error("*** MLS setup failed!\n");
	     }
	   }
	 }
	 else if (pre_or_post == ML_BOTH) {

	   ml->post_smoother[i].data_destroy = ML_Smoother_Destroy_MLS;
	   sprintf(str,"MLS_pre%d",i);

	   ML_Smoother_Set(&(ml->pre_smoother[i]), 
                        (void *) widget, fun, ntimes,
			0.0, str);
	   sprintf(str,"MLS_post%d",i);
	   errCode = ML_Smoother_Set(&(ml->post_smoother[i]), 
               (void *) widget, fun, ntimes, 0.0, str);


	   /* set up the values needed for MLS  */

	   if (fun == ML_Smoother_MLS_Apply) {
	     if (ML_MLS_Setup_Coef(&(ml->post_smoother[i]), 1,ml->symmetrize_matrix)) {
	       return pr_error("*** MLS setup failed!\n");
	     }
	   }
	 }
	 else return(pr_error("Print unknown pre_or_post choice\n"));

#ifdef ML_TIMING
         if (pre_or_post == ML_PRESMOOTHER)
	   ml->pre_smoother[i].build_time = GetClock() - t0;
	 else ml->post_smoother[i].build_time = GetClock() - t0;

         ml->timing->total_build_time   += ml->post_smoother[i].build_time;
#endif
     }
   }
   return errCode;
}
#endif/*MB_MODIF*/

/********************************************************************
1-step krylov method for equivalent real forms. This smoother will
take 1 step of a krylov method applied to D^-1 A where A is in 
equivalent real form and D is the equivalent real form of the diagonal 
complex matrix system.
********************************************************************/
int ML_Gen_Smoother_ERF_1StepKrylov(ML *ml, int nl, int pre_or_post)
{
   int              start_level, end_level, i, errCode=0;
   ML_Operator     *Amat;
   char             str[80];
   int                (*fun)(ML_Smoother *, int, double *, int, double *);
   int ntimes = 1;

#ifdef ML_TIMING
   double         t0;
#endif

   if (nl == ML_ALL_LEVELS) {
	  start_level = 0; end_level = ml->ML_num_levels-1;}
   else { start_level = nl; end_level = nl;}

   if (start_level < 0) {
      printf("ML_Gen_Smoother_ERF_1StepKrylov: cannot set smoother on level %d\n",
		      start_level);
      return 1;
   }

   fun = ML_DiagScaled_1stepKrylov;

   for (i = start_level; i <= end_level; i++) {
#ifdef ML_TIMING
     t0 = GetClock();
#endif
     Amat = &(ml->Amat[i]);
     if (Amat->matvec->func_ptr != NULL) {
	 if (pre_or_post == ML_PRESMOOTHER) {
           sprintf(str,"ERF_1STEP_pre%d",i);
           errCode=ML_Smoother_Set(&(ml->pre_smoother[i]), 
				  NULL, fun, ntimes, 0.0, str);
	 }
	 else if (pre_or_post == ML_POSTSMOOTHER) {
           sprintf(str,"ERF_1STEP_post%d",i);
	   errCode=ML_Smoother_Set(&(ml->post_smoother[i]), 
				  NULL, fun, ntimes, 0.0, str);
	 }
	 else if (pre_or_post == ML_BOTH) {
           sprintf(str,"ERF_1STEP_pre%d",i);

	   ML_Smoother_Set(&(ml->pre_smoother[i]), 
			   NULL, fun, ntimes, 0.0, str);
           sprintf(str,"ERF_1STEP_post%d",i);
	   errCode = ML_Smoother_Set(&(ml->post_smoother[i]), 
				     NULL, fun, ntimes, 0.0, str);
	 }
	 else return(pr_error("Print unknown pre_or_post choice\n"));

#ifdef ML_TIMING
         if (pre_or_post == ML_PRESMOOTHER)
	   ml->pre_smoother[i].build_time = GetClock() - t0;
	 else ml->post_smoother[i].build_time = GetClock() - t0;

         ml->timing->total_build_time   += ml->post_smoother[i].build_time;
#endif
     }
   }
   return errCode;
}



/* ------------------------------------------------------------------------- */
/* generate the sparse approximate inverse smoother */
/* ------------------------------------------------------------------------- */

#ifdef HAVE_ML_PARASAILS
#include "Matrix.h"
#include "ParaSails.h"
#endif

int ML_Gen_Smoother_ParaSails(ML *ml, int nl, int pre_or_post, int ntimes,
   int sym, double thresh, int num_levels, double filter, int parasails_loadbal,
   int parasails_factorized)
{
#ifdef HAVE_ML_PARASAILS
   int            (*fun1)(ML_Smoother *, int, double *, int, double *);
   int            (*fun2)(ML_Smoother *, int, double *, int, double *);
   int            start_level, end_level, i, status = 1;
   int            row, start_row, end_row, row_length;

   Matrix *mat;
   struct widget { int parasails_factorized; ParaSails *ps;} *widget;
   ParaSails *ps;
   int j;
   char str[80];
   
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
	 if( ML_Get_PrintLevel() > 8 ) ParaSailsStatsPattern(ps, mat);
	 ParaSailsSetupValues(ps, mat, filter);
	 if( ML_Get_PrintLevel() > 8 ) ParaSailsStatsValues(ps, mat);

	 /* we can destroy the matrix now */
	 MatrixDestroy(mat);

	 ml->post_smoother[i].data_destroy = ML_Smoother_Clean_ParaSails;

/* Turn on temporarily */
         widget->parasails_factorized = parasails_factorized;
         widget->ps                   = ps;

	 /* by default we want pre- and post-smoother */
	 sprintf(str,"ParaSails_post%d",i);
	 status = ML_Smoother_Set(&(ml->post_smoother[i]), 
				  (void *) widget, fun1, ntimes, 0.0,str);
	 sprintf(str,"ParaSails_pre%d",i);
	 status = ML_Smoother_Set(&(ml->pre_smoother[i]), 
				  (void *) widget, fun2, ntimes, 0.0,str);

#ifdef ML_TIMING
         ml->post_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->post_smoother[i].build_time;
#endif
      }

   /* note: in free, post and pre are the same */

   return(status);
#else
   printf("ParaSails not linked\n");
   ML_avoid_unused_param((void *) ml);
   ML_avoid_unused_param((void *) &nl);
   ML_avoid_unused_param((void *) &pre_or_post);
   ML_avoid_unused_param((void *) &ntimes);
   ML_avoid_unused_param((void *) &sym);
   ML_avoid_unused_param((void *) &thresh);
   ML_avoid_unused_param((void *) &num_levels);
   ML_avoid_unused_param((void *) &filter);
   ML_avoid_unused_param((void *) &parasails_loadbal);
   ML_avoid_unused_param((void *) &parasails_factorized);
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
/*
printf("    (%d): Amat(%d,%d)  Rmat(%d,%d)  Pmat(%d,%d)\n",
                Amat->comm->ML_mypid,
                Amat->outvec_leng,Amat->invec_leng,
                Rmat->outvec_leng,Rmat->invec_leng,
                Pmat->outvec_leng,Pmat->invec_leng);
fflush(stdout);
*/

   if (Amat->matvec->func_ptr == NULL) {
      if (output_level > 3)
      printf("Warning: No Amat matvec on grid %d (where finest = 0).\n\
		can not check Amat's getrow\n",i);
   }
/*
   else ML_Operator_Check_Getrow( Amat,ml->comm,i,"Amat");
*/

   if (Amat->getrow->func_ptr == NULL)
      pr_error("Error: No A matrix getrow on grid %d : \
                       can not do ML_Gen_Amatrix_RAP.\n",i);

   if ((Amat->getrow->pre_comm  == NULL) &&
       (Amat->getrow->post_comm == NULL) && (ml->comm->ML_nprocs > 1) ) {
       if ((ml->comm->ML_mypid == 0) && (output_level > 3)) {
         printf("Warning:No communication information given with Amat's \n");
         printf("\tgetrow on level %d (finest = 0).!!!!\n",i);
       }
   }

   if (Rmat->matvec->func_ptr == NULL) {
      if (output_level > 3)
      printf("Warning: No Rmat matvec on grid %d (where finest = 0).\n\
		can not check Rmat's getrow\n",i);
   }
/*
   else ML_Operator_Check_Getrow( Rmat,ml->comm,i,"Rmat");
*/

   if (Rmat->getrow->func_ptr == NULL)
      pr_error("Error: No R matrix getrow on grid %d : \n\
                       can not do ML_Gen_AmatrixRAP.\n",i);

   if ((Rmat->getrow->pre_comm  == NULL) &&
       (Rmat->getrow->post_comm == NULL) && (ml->comm->ML_nprocs > 1) ) {
       if ((ml->comm->ML_mypid == 0) && (output_level > 3)) {
         printf("Warning:No communication information given with Rmat's \n");
         printf("\tgetrow on level %d (finest = 0).!!!!\n",i);
       }
   }

   if (Pmat->matvec->func_ptr == NULL) {
      if (output_level > 3)
      printf("Warning: No Pmat matvec on grid %d (where finest = 0).\n\
		can not check Pmat's getrow\n",i);
   }
/*
   else ML_Operator_Check_Getrow(Pmat,ml->comm,i,"Pmat");
*/

   if (Pmat->getrow->func_ptr == NULL)
      pr_error("Error: No P matrix getrow on grid %d : \n\
                       can not do ML_Gen_AmatrixRAP.\n",i);

   if ((Pmat->getrow->pre_comm  == NULL) &&
       (Pmat->getrow->post_comm == NULL) && (ml->comm->ML_nprocs > 1) ) {
       if ((ml->comm->ML_mypid == 0) && (output_level > 3)) {
         printf("Warning:No communication information given with Pmat's \n");
         printf("\tgetrow on level %d (finest = 0).!!!!\n",i);
       }
   }

/*ML_Operator_Print(&(ml->Pmat[child_level]),"Pn");*/
   ML_rap(&(ml->Rmat[parent_level]), &(ml->Amat[parent_level]),
          &(ml->Pmat[child_level]), &(ml->Amat[child_level]),
          ML_MSR_MATRIX);

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
                            int (*func)(void*,double*,double*) )
{
   ML_Mapper_SetFunc( ml->SingleLevel[level].eqn2grid,fleng,tleng,func);
   ML_Mapper_SetData( ml->SingleLevel[level].eqn2grid,
                      (void*)ml->SingleLevel[level].Grid->Grid);
   return 0;
}

int ML_Set_GridToEqnMapFunc(ML *ml, int level, int fleng, int tleng,
                            int (*func)(void*,double*,double*) )
{
   ML_Mapper_SetFunc( ml->SingleLevel[level].grid2eqn,fleng,tleng,func);
   ML_Mapper_SetData( ml->SingleLevel[level].grid2eqn,
                      (void*)ml->SingleLevel[level].Grid->Grid);
   return 0;
}

/* ************************************************************************* */
/* Setting up a given solver                                                 */
/* ------------------------------------------------------------------------- */

int ML_Setup(ML *ml, int scheme, int finest_level, int incr_or_decr, void *obj)
{
   ML_Aggregate *aggr;
   ML_AMG       *amg;

   if ( scheme == ML_MG2CGC )
   {
      if ( incr_or_decr == ML_INCREASING )
         ML_Gen_Solver( ml, scheme, finest_level, finest_level+1 );
      else
         ML_Gen_Solver( ml, scheme, finest_level, finest_level-1 );
   }
   else if ( scheme == ML_RSAMG )
   {
      if ( obj != NULL )
      {
         amg = (ML_AMG *) obj;
         if ( amg->ML_id != ML_ID_AMG )
         {
            printf("ML_Setup ERROR : method = RSAMG, data not ML_AMG.\n");
            exit(1);
         }
      }
#ifdef ML_CPP
      ML_Gen_MGHierarchy_UsingAMG(ml, finest_level, incr_or_decr, (ML_AMG *)obj);
#else
      ML_Gen_MGHierarchy_UsingAMG(ml, finest_level, incr_or_decr, obj);
#endif
   }
   else if ( scheme == ML_SAAMG )
   {
      if ( obj != NULL )
      {
         aggr = (ML_Aggregate *) obj;
         if ( aggr->ML_id != ML_ID_AGGRE )
         {
            printf("ML_Setup ERROR : method = SAAMG, data not ML_Aggregate.\n");
            exit(1);
         }
      }
#ifdef ML_CPP
      ML_Gen_MGHierarchy_UsingAggregation(ml, finest_level, incr_or_decr, (ML_Aggregate *)obj);
#else
      ML_Gen_MGHierarchy_UsingAggregation(ml, finest_level, incr_or_decr, obj);
#endif
   }
   else
   {
      printf("ML_Setup ERROR : method not recognized.\n");
      exit(1);
   }
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

   ML_memory_check("gen solver");

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
      if (current_level->Amat->matvec->func_ptr == NULL &&
          level != coarsest_level) {
         pr_error("Error: No A matrix on grid %d.\n",level);
      }

      if ((current_level->Amat->getrow->pre_comm  == NULL) &&
          (current_level->Amat->getrow->post_comm == NULL) &&
          (current_level->Amat->getrow->func_ptr != NULL) &&
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
         if (current_level->Rmat->matvec->func_ptr == NULL)
            pr_error("Error: No R matvec on grid %d.\n",level);
/*
         ML_Operator_Check_Getrow(current_level->Rmat,ml->comm,level,"Rmat");
*/
         if (level != finest_level &&
             current_level->Pmat->matvec->func_ptr == NULL)
            pr_error("Error: No P matvec on grid %d.\n",level);
/*
         ML_Operator_Check_Getrow(current_level->Pmat,ml->comm,level, "Pmat");
*/
      }

      if (current_level->pre_smoother->smoother->func_ptr==ML_Smoother_Jacobi)
      {
         if ((temp == NULL) && (current_level->csolve->func->func_ptr == NULL))
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
         ML_gsum_scalar_int(&t1, &j, ml->comm);
         ML_gsum_scalar_int(&t2, &j, ml->comm);
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
         dtmp1 = (double*) ML_allocate( leng1 * sizeof(double) );
         dtmp2 = (double*) ML_allocate( leng2 * sizeof(double) );
         ML_BdryPts_Get_Dirichlet_Grid_Info(current_level->BCs,&leng3,&itmp3);
         for ( j = 0; j < leng1; j++ ) dtmp1[j] = 0.0;
         for ( j = 0; j < leng2; j++ ) dtmp2[j] = 0.0;
         for ( j = 0; j < leng3; j++ ) dtmp1[itmp3[j]] = 1.0;
         ML_Mapper_Apply(current_level->grid2eqn, dtmp1, dtmp2);
         leng1 = 0;
         for ( j = 0; j < leng2; j++ ) if ( dtmp2[j] == 1.0 ) leng1++;
         itmp3 = (int*) ML_allocate( leng1 * sizeof(int) );
         leng1 = 0;
         for ( j = 0; j < leng2; j++ )
            if ( dtmp2[j] == 1.0 ) itmp3[leng1++] = j;
         ML_BdryPts_Load_Dirichlet_Eqn(current_level->BCs, leng1, itmp3);
         ML_free( itmp3 );
      } else {
         ML_BdryPts_Copy_Dirichlet_GridToEqn(current_level->BCs);
      }
      current_level = temp;
   }
   ml->ML_num_actual_levels = i;
   /*
   if ((ml->comm->ML_mypid == 0) && (output_level > 0))
      printf("Total number of levels = %d\n",i);
   */

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
/* ML_Solve : call iteration function                                        */
/*-------------------------------------------------------------------------- */

int ML_Solve(ML *ml, int inlen, double *sol, int outlen, double *rhs)
{
   int level, leng;

   level = ml->ML_finest_level;
   leng  = ml->Amat[level].outvec_leng;
   if ( inlen != outlen || leng != inlen )
   {
      printf("ML_Solve ERROR : inlen != outlen != mat length.\n");
      exit(1);
   }
   ML_Iterate(ml, sol, rhs);
   return 0;
}

/* ************************************************************************* */
/* ML_Iterate : call iteration function                                      */
/*-------------------------------------------------------------------------- */

int ML_Iterate(ML *ml, double *sol, double *rhs)
{
   int  i = 0, count = 0;
   double res_norm, prev_res_norm = -1.0, reduction, r0 = 0.0
	  /*, old_reduction = 10.*/;

   reduction = 1.;

   if ((ml->res_output_freq > 0) && (ml->comm->ML_mypid == 0) ) {
      printf("Iter       ||res_i||_2**    ||res_i||/||res_i+1||\n");
      count = ml->res_output_freq - 1;
   }
   else count = 0;

   while ( reduction >= ml->tolerance && i < ml->max_iterations )
   {
      res_norm = ML_Cycle_MG( &(ml->SingleLevel[ml->ML_finest_level]), sol,
				rhs, ML_NONZERO,ml->comm, ML_COMPUTE_RES_NORM,ml);
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
      if (i == 1) r0 = res_norm + 1.e-25;
      prev_res_norm = res_norm + 1.e-25;
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
/* set to 1 to use another ML cycle after ggb, to 0 otherwise                */
/*-------------------------------------------------------------------------- */

static int ML_ggb_SymmetricCycle = 1;
int ML_ggb_Set_SymmetricCycle(int flag) 
{
  ML_ggb_SymmetricCycle = flag;
  return 0;
}

/*****************************************************************************/
/* set to 1 to use First GGB cycle and then ML cycle, to 0 for ML first      */
/*-------------------------------------------------------------------------- */
static int ML_GGBcycFirst = 0;
static int store_AQ = 0;       /* May not be defined if storage is an issue. Works only with GGBcycFirst
                                  and SuperLU   */
int ML_ggb_Set_GGBCycle(int flag) 
{
  ML_GGBcycFirst = flag;
  store_AQ = flag;
  return 0;
}



/****************************************************/
/* Switch between SuperLU (default) and Amesos      */
/* Note: store AQ option works only with SuperLU    */
/*       (For Now...) and   ML_GGBcycFirst = 1      */ 
/*------------------------------------------------- */
#if defined(HAVE_ML_ARPACK) || defined(HAVE_ML_PARPACK)
#if defined(HAVE_ML_SUPERLU) || defined(HAVE_ML_AMESOS)
#ifdef HAVE_ML_SUPERLU
static int ML_ggb_CoarseSolver = 1;
#else
static int ML_ggb_CoarseSolver = 2;
#endif

#else
 fprintf( stderr,
	   "ERROR: ML has not been configured with either AMESOS or SUPERLU support.\n"
	   "ERROR: Please reconfigure.\n"
	   "ERROR: (file %s, line %d)\n",
	   __FILE__, __LINE__ );
  exit(-1);
#endif

#else
static int ML_ggb_CoarseSolver = 1;
#endif /* if defined(HAVE_ML_ARPACK) || defined(HAVE_ML_PARPACK) */

int ML_ggb_Set_CoarseSolver(int flag) 
{
  ML_ggb_CoarseSolver = flag;
  return 0;
}

/*****************************************************************************/
/* solve using V-cycle multigrid                                             */
/*-------------------------------------------------------------------------- */
int ML_Solve_MGV( ML *ml , double *din, double *dout)
{
  int    i, leng, dir_leng, *dir_list, k, level;
  double *diag, *scales, *din_temp, *dout_tmp;
  ML     *ml_ggb;
 
   /* ------------------------------------------------------------ */
   /* initially set the solution to be all 0           	           */
   /* ------------------------------------------------------------ */

   level = ml->ML_finest_level;
   leng = ml->Amat[level].outvec_leng;
   for ( i = 0; i < leng; i++ ) dout[i] = 0.0;
   din_temp = (double*) ML_allocate( leng * sizeof(double) );

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

   if (ml->void_options == NULL)
     ML_Cycle_MG(&(ml->SingleLevel[ml->ML_finest_level]), dout, din_temp,
		 ML_ZERO, ml->comm, ML_NO_RES_NORM, ml);
   
   else
     {
       
       ml_ggb = (ML *) ml->void_options;
       

       
       /* Notify user */
       if (ML_GGBcycFirst == 1 && ML_ggb_SymmetricCycle == 1) {
	 ML_ggb_SymmetricCycle =0;
	 printf(" Symmetric cycle goes only when GGB is the second cycle \n");
	 printf(" Switching to Nonsymmetric mode     \n"); 	 
       }


       if (ML_GGBcycFirst == 1) { /* First Do GGB cycle then MG */ 

	 /* Option for SuperLU solver */
	 if (ML_ggb_CoarseSolver == 1) {
	   double    *sol;
	   
	   sol  = (double *) ML_allocate(leng*sizeof(double));
	   
	   /* Manualy perform the 2 level GGB cycle */       
	   ML_Cycle_GGB(ml_ggb, sol, din_temp);
	   
	   ML_Cycle_MG(&(ml->SingleLevel[ml->ML_finest_level]), dout, din_temp,
		       ML_ZERO, ml->comm, ML_NO_RES_NORM, ml);
	   
	   /* Add both solutions to get the right one */ 
	   for ( i = 0; i < leng; i++ ) dout[i] += sol[i];    	 
	   
	   ML_free(sol);
	   
	 }
	 else /* Option for Amesos solver */
	   {
	     ML_Cycle_MG(&(ml_ggb->SingleLevel[ml_ggb->ML_finest_level]), dout,
			 din_temp, ML_ZERO, ml_ggb->comm, ML_NO_RES_NORM, ml_ggb);	   
	     
	     ML_Cycle_MG(&(ml->SingleLevel[ml->ML_finest_level]), dout, din_temp,
			 ML_NONZERO, ml->comm, ML_NO_RES_NORM, ml);
	     
	   }		       
       } 
       else  /* Do GGB cycle second */
	 {
	 /* First do the MG part */
	 ML_Cycle_MG(&(ml->SingleLevel[ml->ML_finest_level]), dout, din_temp,
		     ML_ZERO, ml->comm, ML_NO_RES_NORM, ml);
	 
	 
	 if (ML_ggb_CoarseSolver == 1) {
	   /* Manualy perform the 2 level GGB cycle */
	   ML_Cycle_GGB(ml_ggb, dout, din_temp);
	 }
	 else  /* Option for Amesos solver */
	   {
	     ML_Cycle_MG(&(ml_ggb->SingleLevel[ml_ggb->ML_finest_level]), dout,
			 din_temp, ML_NONZERO, ml_ggb->comm, ML_NO_RES_NORM, ml_ggb);	 
	   }
	 
	 
	 /* Used for symmetric GGB cycle */
	 
	 if (ML_ggb_SymmetricCycle == 1)
	   ML_Cycle_MG(&(ml->SingleLevel[ml->ML_finest_level]), dout, din_temp,
		       ML_NONZERO, ml->comm, ML_NO_RES_NORM, ml);
	 
	 
	 
	 }
     }
   
       
       /*	 
       ML_Cycle_MG(&(ml_ggb->SingleLevel[ml_ggb->ML_finest_level]), dout,
       din_temp, ML_NONZERO, ml_ggb->comm, ML_NO_RES_NORM, ml_ggb);
       */          
              
       /* Trying the additive version (== doesn't work) */
       /*     
	      dout_tmp = (double*) ML_allocate( leng * sizeof(double) );
	      for ( i = 0; i < leng; i++ )        dout_tmp[i] = 0.0;
	      
	      ML_Cycle_MG(&(ml_ggb->SingleLevel[ml_ggb->ML_finest_level]), dout_tmp,
	      din_temp, ML_ZERO, ml_ggb->comm, ML_NO_RES_NORM, ml_ggb);
	      
	      for ( i = 0; i < leng; i++ )  dout[i] = dout[i] + dout_tmp[i];
	      ML_free(dout_tmp);
       */
       

       
       /* In the cheap cycle we do the two level cycle manualy saving a mat vec */
       
       
       /* "after cycle" only if required (default is yes) */



   
   
   /*
   ML_Cycle_MG(&(ml->SingleLevel[ml->ML_finest_level]), dout, din_temp,
	       ML_NONZERO, ml->comm, ML_NO_RES_NORM, ml);
   */
   ML_free(din_temp);
   return 0;
}

/*****************************************************************************/
/* 2 Level GGB cycle performed manualy to save mat-vec                       */
/*-------------------------------------------------------------------------- */

extern int ML_Cycle_GGB(ML *ml_ggb, double *sol, double *rhs)
{
  

  int      lengc, lengf, i;
  double  *rhs0, *sol0; 
  ML_Operator *Rmat, *Qtilde, *Pmat;
  

  
  Rmat   = &(ml_ggb->Rmat[1]);
  Pmat   = &(ml_ggb->Pmat[0]);
  
  
  /* Haim Trying Dense Mat vec */
  lengc = Rmat->outvec_leng;
  lengf = Pmat->outvec_leng;
  
  rhs0  = (double *) ML_allocate(lengc*sizeof(double));
  sol0  = (double *) ML_allocate(lengc*sizeof(double));
  

  if (ML_GGBcycFirst == 1) {  /* Do GGB cycle First and later MG */

    double  *tmp1;
    tmp1  = (double *) ML_allocate(lengf*sizeof(double));
    
    /* get Q'*b */
    ML_Operator_Apply(Rmat, lengf, rhs, lengc, rhs0);
    
    /* Solve coarse grid problem */
    ML_CSolve_Apply( Rmat->to->csolve, lengc, sol0, lengc, rhs0);  
    
    
    /* Need to do it anyway */
    ML_Operator_Apply(Pmat, lengc, sol0, lengf, sol);
      
    
    /* Use this option if A*Q was stored in the process of computing rap in GGB build */    
    if (store_AQ == 1) {
      
      /* Qtilde = Q'*A  qas stored in the coarse grid build */
      Qtilde = (ML_Operator *) ml_ggb->void_options; 
      
      /* Interpolate to fine mesh */ /* Here we should put Qtilde */
      ML_Operator_Apply(Qtilde, lengc, sol0, lengf, tmp1);
    }
    else 
      
      /* Compute K*u . Can drop this line after we put Qtilde*/ 
      ML_Operator_Apply( &(ml_ggb->Amat[1]), lengf, sol, lengf, tmp1);
    
    
    /* Correct the residual term */
    for ( i = 0; i < lengf; i++ ) rhs[i] = rhs[i] - tmp1[i];
    
    ML_free(tmp1);
    
  }
  
  else /* If GGB cycle is second */
    {

      /* # ifdef GGBcycSecond */
       /* We need to compute the following coarse grid correction:  
        *	  u = u + Q*(Ac)^(-1) *Q'A*r =
	*           = u + Q*(Ac)^(-1) *(Qtilde*u - Q'*b)	  
       */  
      double *tmp1, *tmp2, *rhs_tmp;
      tmp1  = (double *) ML_allocate(lengc*sizeof(double));
      tmp2  = (double *) ML_allocate(lengf*sizeof(double));
      rhs_tmp = (double *) ML_allocate(lengf*sizeof(double));
      
#ifdef store_QtransA       
      
      /* Qtilde = Q'*A  qas stored in the coarse grid build */
      Qtilde = (ML_Operator *) ml_ggb->void_options; 
      
      
      /* get Q'*b */
      ML_Operator_Apply(Rmat, lengf, rhs, lengc, rhs0);
      
      
      /* get Qtilde*u */
      ML_Operator_Apply(Qtilde, lengf, sol, lengc, tmp1);
      
      /* compute   Qtilde*u - Q'*b. which forms the right hand side of the coarse grid */
      for ( i = 0; i < lengc; i++ )  rhs0[i] =  rhs0[i] - tmp1[i];
      
#else 
      
      /* first compute K*u */
      ML_Operator_Apply( &(ml_ggb->Amat[1]), lengf, sol, lengf, tmp2);
      
      /* Get the fine grid residual */
      for ( i = 0; i < lengf; i++ )  rhs_tmp[i] =  rhs[i] - tmp2[i];
      
      /* Compute the rhs on the coarse grid */
      ML_Operator_Apply(Rmat, lengf, rhs_tmp, lengc, rhs0);
      
#endif  
      
      /* Solve coarse grid problem */
      ML_CSolve_Apply( Rmat->to->csolve, lengc, sol0, lengc, rhs0);  
      
      /* Interpolate to fine mesh */
      ML_Operator_Apply(Pmat, lengc, sol0, lengf, tmp2);
      
      /* Correct the solution */
      for ( i = 0; i < lengf; i++ ) sol[i] += tmp2[i];
      
      
      ML_free(tmp1);
      ML_free(tmp2);
      ML_free(rhs_tmp);
      
    }







  ML_free(sol0);
  ML_free(rhs0);
 
  return 1;
}

/*****************************************************************************/
/* solve using Full MGV-cycle multigrid                                             */
/*-------------------------------------------------------------------------- */

int ML_Solve_MGFull( ML *ml , double *din, double *dout)
{
   int    i, leng, dir_leng, *dir_list, k, level;
   double *diag, *scales, *din_temp;

   /* ------------------------------------------------------------ */
   /* initially set the solution to be all 0           	           */
   /* ------------------------------------------------------------ */

   level = ml->ML_finest_level;
   leng = ml->Amat[level].outvec_leng;
   for ( i = 0; i < leng; i++ ) dout[i] = 0.0;
   din_temp = (double*) ML_allocate( leng * sizeof(double) );

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

   ML_Cycle_MGFull(&(ml->SingleLevel[ml->ML_finest_level]), dout, din_temp,
                ML_ZERO, ml->comm, ML_NO_RES_NORM,ml);

   ML_free(din_temp);
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
   if (res == NULL) {
     pr_error("swillie: out of space\n");
     ML_avoid_unused_param((void *) &isize);
     ML_avoid_unused_param((void *) &osize);
   }

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
   din_temp = (double*) ML_allocate( leng * sizeof(double) );

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

   ML_Cycle_MG(&(ml->SingleLevel[ml->ML_finest_level]), dout, din_temp,
                ML_ZERO, ml->comm, ML_NO_RES_NORM, ml);

   ML_free(din_temp);
   return 0;
}

#ifdef ML_SINGLE_LEVEL_PROFILING

/* prototypes for stubs */
double ML_Cycle_MG_9(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml);
double ML_Cycle_MG_8(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml);
double ML_Cycle_MG_7(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml);
double ML_Cycle_MG_6(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml);
double ML_Cycle_MG_5(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml);
double ML_Cycle_MG_4(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml);
double ML_Cycle_MG_3(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml);
double ML_Cycle_MG_2(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml);
double ML_Cycle_MG_1(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml);
double ML_Cycle_MG_0(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml);
/* end of stub prototypes */

/* Stub functions.  The include file "ml_v_cycle.c" has the actual source
   code.  This is the same code as found in function ML_Cycle_MG() below. */
double ML_Cycle_MG_9(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
#include "ml_v_cycle.c"
}

double ML_Cycle_MG_8(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
#include "ml_v_cycle.c"
}

double ML_Cycle_MG_7(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
#include "ml_v_cycle.c"
}

double ML_Cycle_MG_6(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
#include "ml_v_cycle.c"
}

double ML_Cycle_MG_5(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
#include "ml_v_cycle.c"
}

double ML_Cycle_MG_4(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
#include "ml_v_cycle.c"
}

double ML_Cycle_MG_3(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
#include "ml_v_cycle.c"
}

double ML_Cycle_MG_2(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
#include "ml_v_cycle.c"
}

double ML_Cycle_MG_1(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
#include "ml_v_cycle.c"
}

double ML_Cycle_MG_0(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
#include "ml_v_cycle.c"
}
#endif    /* ifdef ML_SINGLE_LEVEL_PROFILING */

/*****************************************************************************/
/* solve using V-cycle multigrid                                             */
/*-------------------------------------------------------------------------- */

double ML_Cycle_MG(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
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

   rhss = (double *) ML_allocate( lengf * sizeof(double) );
   ML_DVector_GetDataPtr(curr->Amat_Normalization, &normalscales) ;
   for ( i = 0; i < lengf; i++ ) rhss[i] = rhs[i];

#ifdef ML_ANALYSIS
   if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
   {
      fp = fopen("mlmatlab.m", "w");
      Nrows = lengf;
      fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
      allocated_space = 100;
      cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
      vals = (double *) ML_allocate(allocated_space*sizeof(double));
      for (i = 0; i < lengf; i++) {
         while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)==0)
         {
            allocated_space = 2*allocated_space + 1;
            ML_free(vals); ML_free(cols);
            cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
            vals = (double *) ML_allocate(allocated_space*sizeof(double));
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
         res = (double *) ML_allocate(lengf*sizeof(double));
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhss[i] - res[i];
         res_norm = sqrt(ML_gdot(lengf, res, res, comm));
         ML_free(res);
      }
   }
   else {
      res = (double *) ML_allocate(lengf*sizeof(double));

      /* --------------------------------------------------------- */
      /* pre-smoothing and compute residual                        */
      /* --------------------------------------------------------- */
      ML_Smoother_Apply(pre, lengf, sol, lengf, rhss, approx_all_zeros);

      if ( ( approx_all_zeros != ML_ZERO ) ||
           ( pre->smoother->func_ptr != NULL) )
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
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
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
         squareA = ML_allocate( Nrows * Nrows * sizeof(double) );
         for ( i = 0; i < Nrows*Nrows; i++ ) squareA[i] = 0.0;
         allocated_space = 100;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               squareA[cols[j]*Nrows+i] = vals[j];
         }
         ML_free(cols); ML_free(vals);
         eig = (double *) ML_allocate( Nrows * sizeof(double) );
         work = (double *) ML_allocate( 4 * Nrows * sizeof(double) );
         lwork = 4 * Nrows;
         jobz = 'V'; jobz2 = 'U';
         MLFORTRAN(dsyev)(&jobz,&jobz2,&Nrows,squareA,&Nrows,eig,work,&lwork,&info,
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
         ML_free(squareA);
         ML_free(eig);
         ML_free(work);
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

      rhs2 = (double *) ML_allocate(lengc*sizeof(double));
      sol2 = (double *) ML_allocate(lengc*sizeof(double));
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
         dtmp = (double *) ML_allocate( lengf * sizeof( double ) );
         ML_Mapper_Apply(curr->eqn2grid, res, dtmp );
         for ( i = 0; i < lengf; i++ ) res[i] = dtmp[i];
         ML_free( dtmp );
      }
      ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, res, lengc, rhs2);
      if ( ML_Mapper_Check(Rmat->to->grid2eqn) == 1 )
      {
         dtmp = (double *) ML_allocate( lengc * sizeof( double ) );
         ML_Mapper_Apply(Rmat->to->grid2eqn, rhs2, dtmp );
         for ( i = 0; i < lengc; i++ ) rhs2[i] = dtmp[i];
         ML_free( dtmp );
      }
      ML_DVector_GetDataPtr(Rmat->to->Amat_Normalization,&normalscales);
      if ( normalscales != NULL )
         for ( i = 0; i < lengc; i++ ) rhs2[i] = rhs2[i] * normalscales[i];

      /* --------------------------------------------------------- */
      /* process the next level and transfer back to this level    */
      /* --------------------------------------------------------- */
#ifdef ML_SINGLE_LEVEL_PROFILING
      i = curr->levelnum-1;
      switch(i) {
         case 9:
      ML_Cycle_MG_9( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_9( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 8:
      ML_Cycle_MG_8( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_8( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 7:
      ML_Cycle_MG_7( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_7( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 6:
      ML_Cycle_MG_6( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_6( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 5:
      ML_Cycle_MG_5( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_5( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 4:
      ML_Cycle_MG_4( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_4( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 3:
      ML_Cycle_MG_3( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_3( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 2:
      ML_Cycle_MG_2( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_2( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 1:
      ML_Cycle_MG_1( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_1( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 0:
      ML_Cycle_MG_0( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_0( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         default:
      ML_Cycle_MG( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
   }

#else
      ML_Cycle_MG( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
#endif /* ifdef ML_SINGLE_LEVEL_PROFILING */

      /* ------------------------------------------------------------ */
      /* transform the data from equation to grid space, do grid      */
      /* transfer and then transfer back to equation space            */
      /* ------------------------------------------------------------ */
      if ( ML_Mapper_Check(Rmat->to->eqn2grid) == 1 )
      {
         dtmp = (double *) ML_allocate( lengc * sizeof( double ) );
         ML_Mapper_Apply(Rmat->to->eqn2grid, sol2, dtmp);
         for ( i = 0; i < lengc; i++ ) sol2[i] = dtmp[i];
         ML_free( dtmp );
      }
      ML_Operator_ApplyAndResetBdryPts(Rmat->to->Pmat,lengc,sol2,lengf,res);
      if ( ML_Mapper_Check(curr->grid2eqn) == 1 )
      {
         dtmp = (double *) ML_allocate( lengf * sizeof( double ) );
         ML_Mapper_Apply(curr->grid2eqn, res, dtmp);
         for ( i = 0; i < lengf; i++ ) res[i] = dtmp[i];
         ML_free( dtmp );
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
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
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
      ML_free(sol2);
      ML_free(rhs2);

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
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
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
      ML_free(res);
   }

   ML_free(rhss);
   return(res_norm);
}
/*****************************************************************************/
/* solve using V-cycle multigrid                                             */
/*-------------------------------------------------------------------------- */

double ML_Cycle_MGFull(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{
   int         i, lengc, lengf;
   double      *res,  *sol2, *rhs2, res_norm = 0., *normalscales;
   double      *rhss, *dtmp;
   ML_Operator *Amat, *Rmat;

   Amat     = curr->Amat;
   Rmat     = curr->Rmat;
   lengf    = Amat->outvec_leng;

   /* ------------------------------------------------------------ */
   /* first do the normalization                                   */
   /* ------------------------------------------------------------ */

   rhss = (double *) ML_allocate( lengf * sizeof(double) );
   ML_DVector_GetDataPtr(curr->Amat_Normalization, &normalscales) ;
   for ( i = 0; i < lengf; i++ ) rhss[i] = rhs[i];

   /* ------------------------------------------------------------ */
   /* smoothing or coarse solve                                    */
   /* ------------------------------------------------------------ */
   if (Rmat->to != NULL) {    /* not coarsest grid */
      res = (double *) ML_allocate(lengf*sizeof(double));
      if ( approx_all_zeros != ML_ZERO ) {
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhss[i] - res[i];
      }
      else for ( i = 0; i < lengf; i++ ) res[i] = rhss[i];

      /* project down and do a full cycle */

      lengc = Rmat->outvec_leng;
      rhs2 = (double *) ML_allocate(lengc*sizeof(double));
      sol2 = (double *) ML_allocate(lengc*sizeof(double));
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
         dtmp = (double *) ML_allocate( lengf * sizeof( double ) );
         ML_Mapper_Apply(curr->eqn2grid, res, dtmp );
         for ( i = 0; i < lengf; i++ ) res[i] = dtmp[i];
         ML_free( dtmp );
      }
      ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, res, lengc, rhs2);
      if ( ML_Mapper_Check(Rmat->to->grid2eqn) == 1 )
      {
         dtmp = (double *) ML_allocate( lengc * sizeof( double ) );
         ML_Mapper_Apply(Rmat->to->grid2eqn, rhs2, dtmp );
         for ( i = 0; i < lengc; i++ ) rhs2[i] = dtmp[i];
         ML_free( dtmp );
      }
      ML_DVector_GetDataPtr(Rmat->to->Amat_Normalization,&normalscales);
      if ( normalscales != NULL )
         for ( i = 0; i < lengc; i++ ) rhs2[i] = rhs2[i] * normalscales[i];

      ML_Cycle_MGFull( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);

      /* bring it back up */

      /* ------------------------------------------------------------ */
      /* transform the data from equation to grid space, do grid      */
      /* transfer and then transfer back to equation space            */
      /* ------------------------------------------------------------ */
      if ( ML_Mapper_Check(Rmat->to->eqn2grid) == 1 )
      {
         dtmp = (double *) ML_allocate( lengc * sizeof( double ) );
         ML_Mapper_Apply(Rmat->to->eqn2grid, sol2, dtmp);
         for ( i = 0; i < lengc; i++ ) sol2[i] = dtmp[i];
         ML_free( dtmp );
      }
      ML_Operator_ApplyAndResetBdryPts(Rmat->to->Pmat,lengc,sol2,lengf,res);
      if ( ML_Mapper_Check(curr->grid2eqn) == 1 )
      {
         dtmp = (double *) ML_allocate( lengf * sizeof( double ) );
         ML_Mapper_Apply(curr->grid2eqn, res, dtmp);
         for ( i = 0; i < lengf; i++ ) res[i] = dtmp[i];
         ML_free( dtmp );
      }

      for ( i = 0; i < lengf; i++ ) sol[i] += res[i];
      ML_free(sol2);
      ML_free(rhs2);
      ML_free(res);
      approx_all_zeros = ML_NONZERO;
   }
   ML_free(rhss);
   res_norm = ML_Cycle_MG(curr,sol,rhs,approx_all_zeros,comm,res_norm_or_not, ml);

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

/******************************************************************************
   Solve by first projecting out certain user-defined modes and then calling
   the standard algebraic V-cycle.  We implicitly assume that we are using
   CG as the outer solve, so the projection step is also done after the
   multigrid call to make the entire preconditioner symmetric.

   Note!  We assume that there are no more than 3 modes to project out.
   If there are more, then the constant ML_MAX_SUBSPACE_DIM and loop
   lengths in varios places need to be changed.

   Here is the algorithm, where V' = transpose(V).

        1) Solve V'*A*V * x1 = V' * f.
        2) Calculate a residual r1 = f - A*V*x1.

        3) Solve A * x2 = r1 via a multigrid cycle.

        4) Calculate a residual r2 = f - A * (V*x1 + x2).
        5) Solve V'*A*V * x3 = V'*r2.
        6) Assemble the final solution, x = V*x1 + x2 + V*x3.
******************************************************************************/
int ML_Solve_ProjectedAMGV( ML *ml , double *din, double *dout)
{
   ML_Operator *Amat;
   double **V;
   int dimV, lengV, *pivots;
   double **Av, *VAV;
   int i,j;
   double rhs[ML_MAX_SUBSPACE_DIM], x1[ML_MAX_SUBSPACE_DIM],
          x3[ML_MAX_SUBSPACE_DIM];
   double *res1,*res2;
   double *vec1, *vec2;
   int info;
   /* lapack parameters */
   char trans[2];
   int nrhs=1;
   /*
   unsigned int itmp=0;
   */

   Amat = &(ml->Amat[ml->ML_finest_level]);
   V = Amat->subspace->basis_vectors;
   dimV = Amat->subspace->dimension;
   lengV = Amat->subspace->vecleng;
   VAV   = Amat->subspace->VAV;
   pivots = Amat->subspace->pivots;

   assert(Amat->invec_leng == lengV);

   /* set up V'*A*V if necessary */
   if (Amat->subspace->VAVdone == 0)
   {
      Av = (double **) ML_allocate( dimV * sizeof(double *) );
      for (i=0; i<dimV; i++)
        Av[i] = (double *) ML_allocate( lengV * sizeof(double) );

      for (i=0; i<dimV; i++)
         ML_Operator_Apply(Amat, Amat->invec_leng,V[i],Amat->outvec_leng,Av[i]);
      for (i=0; i<dimV; i++)
        for (j=0; j<dimV; j++) {
          VAV[i*dimV+j] = ML_gdot(lengV, Av[i], V[j], ml->comm);
        }

      /* clean up */
      for (i=0; i<dimV; i++)
        ML_free(Av[i]);
      ML_free(Av);

      /* factor VAV */
      /* see man page for description of arguments */
      DGETRF_F77(&dimV,&dimV,VAV,&dimV,pivots,&info);
      if (info < 0) {
        printf("ML_Solve_ProjectedAMGV: %dth argument to dgetrf has ",-info);
        printf("illegal value\n");
        abort();
      }
      else if (info > 0) {
        printf("ML_Solve_ProjectedAMGV: ");
        printf("U factor is singular in %d'th position\n",i);
        abort();
      }

      Amat->subspace->VAVdone = 1;

   } /* if Amat->subspace->VAVdone == 0 */

   /* calculate f = V^T * b */
   for (i=0; i<dimV; i++)
     rhs[i] = ML_gdot(lengV, V[i], din, ml->comm);

   /**************************************
   **** 1) Solve V'*A*V * x1 = V' * f.
   **************************************/
   /* see man page for description of arguments */
   strcpy(trans,"N");
   /*DGETRS_F77(trans,&dimV,&nrhs,VAV,&dimV,pivots,rhs,&dimV,&info,itmp);*/
   DGETRS_F77(trans,&dimV,&nrhs,VAV,&dimV,pivots,rhs,&dimV,&info);
   if (info < 0) {
     printf("ML_Solve_ProjectedAMGV: %dth argument to dgetrs has ",-info);
     printf("illegal value\n");
     abort();
   }
   for (i=0; i<dimV; i++) x1[i] = rhs[i];

   /*********************************************
   **** 2) Calculate a residual r1 = f - A*V*x1.
   *********************************************/

   res1 = Amat->subspace->res1;
   res2 = Amat->subspace->res2;
   vec1 = Amat->subspace->vec1;
   vec2 = Amat->subspace->vec2;
   /* vec1 = V*x1 */
   for (i=0; i<lengV; i++) {
     vec1[i] = 0.0;
     for (j=0; j<dimV; j++)
        vec1[i] += x1[j] * V[j][i];
   }
   ML_Operator_Apply(Amat, Amat->invec_leng, vec1, Amat->outvec_leng, res1);
   for (i=0; i<Amat->outvec_leng; i++)
     res1[i] = din[i] - res1[i];


   /************************************************
   **** 3) Solve A * x2 = r1 via a multigrid cycle.
   ************************************************/

   ML_Solve_AMGV(ml, res1, dout);

   /******************************************************
   **** 4) Calculate a residual r2 = f - A * (V*x1 + x2).
   **** Here, dout = x2.
   ******************************************************/

   /* vec2 = V*x1 + x2 */
   for (i=0; i<lengV; i++)
     vec2[i] = vec1[i] + dout[i];
   ML_Operator_Apply(Amat, Amat->invec_leng, vec2, Amat->outvec_leng, res2);
   for (i=0; i<lengV; i++)
     res2[i] = din[i] - res2[i];


   /**************************************
   **** 5) Solve V'*A*V * x3 = V'*r2.
   **************************************/

   /* calculate V' * r2 */
   for (i=0; i<dimV; i++)
     rhs[i] = ML_gdot(lengV, V[i], res2, ml->comm);

   /* see man page for description of arguments */
   /*DGETRS_F77(trans,&dimV,&nrhs,VAV,&dimV,pivots,rhs,&dimV,&info,itmp);*/
   DGETRS_F77(trans,&dimV,&nrhs,VAV,&dimV,pivots,rhs,&dimV,&info);
   if (info < 0) {
     printf("ML_Solve_ProjectedAMGV: %dth argument to dgetrs has ",info);
     printf("illegal value\n");
     abort();
   }
   for (i=0; i<dimV; i++) x3[i] = rhs[i];

   /**********************************************************
   **** 6) Assemble the final solution, x = V*x1 + x2 + V*x3.
   **********************************************************/

   /* calculate V * x3 */
   for (i=0; i<lengV; i++) {
     vec1[i] = 0.0;
     for (j=0; j<dimV; j++)
        vec1[i] = vec1[i] + x3[j] * V[j][i];
   }

   /* Recall that vec2 = V*x1 + x2  and vec1 = V*x3 */
   for (i=0; i<Amat->outvec_leng; i++)
     dout[i] = vec2[i] + vec1[i];

   return 0;
}

/*****************************************************************************/
/* function to perform one AMG V-cycle                                       */
/* ------------------------------------------------------------------------- */

double ML_Cycle_AMGV(ML_1Level *curr, double *sol, double *rhs,
	int approx_all_zeros, ML_Comm *comm)
{
   int         i, lengc, lengf;
   double      *res,  *sol2 = NULL, *rhs2 = NULL, res_norm = 0.;
   ML_Operator *Amat, *Rmat;
   ML_Smoother *pre,  *post;
   ML_CSolve   *csolve;
   static int fine_size = 0;

#ifdef RAP2_CHECK
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
      cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
      vals = (double *) ML_allocate(allocated_space*sizeof(double));
      for (i = 0; i < lengf; i++) {
         while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)==0)
         {
            allocated_space = 2*allocated_space + 1;
            ML_free(vals); ML_free(cols);
            cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
            vals = (double *) ML_allocate(allocated_space*sizeof(double));
            if (vals == NULL) {
               printf("Not enough space to get matrix row. Row length of\n");
               printf("%d was not sufficient\n",(allocated_space-1)/2);
               exit(1);
            }
         }
         for (j = 0; j < ncols; j++)
            fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
      }
      ML_free(cols);
      ML_free(vals);
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
         res = (double *) ML_allocate(lengf*sizeof(double));
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhs[i] - res[i];
         res_norm = sqrt(ML_gdot(lengf, res, res, comm));
         ML_free(res);
      }
   }
   else
   {
      res = (double *) ML_allocate(lengf*sizeof(double));

      /* --------------------------------------------------------- */
      /* pre-smoothing and compute residual                        */
      /* --------------------------------------------------------- */

      ML_Smoother_Apply(pre, lengf, sol, lengf, rhs, approx_all_zeros);

      if ( ( approx_all_zeros != ML_ZERO ) ||
           ( pre->smoother->func_ptr != NULL) )
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
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
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
         squareA = ML_allocate( Nrows * Nrows * sizeof(double) );
         for ( i = 0; i < Nrows*Nrows; i++ ) squareA[i] = 0.0;
         allocated_space = 100;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++)
         {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               squareA[cols[j]*Nrows+i] = vals[j];
         }
         ML_free(cols); ML_free(vals);
         eig = (double *) ML_allocate( Nrows * sizeof(double) );
         work = (double *) ML_allocate( 4 * Nrows * sizeof(double) );
         lwork = 4 * Nrows;
         jobz = 'V'; jobz2 = 'U';
         MLFORTRAN(dsyev)(&jobz,&jobz2,&Nrows,squareA,&Nrows,eig,work,&lwork,&info,
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
         ML_free(squareA);
         ML_free(eig);
         ML_free(work);
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
         rhs2 = (double *) ML_allocate(lengc*sizeof(double));
         sol2 = (double *) ML_allocate(lengc*sizeof(double));
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
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
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
      if ( lengc > 0 ) ML_free(sol2);
      if ( lengc > 0 ) ML_free(rhs2);

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
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
         }
         ML_free(cols); ML_free(vals);
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
      ML_free(res);
   }
   return(res_norm);
}

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

     if ((inmat->comminfo != NULL)) {
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
   /*ML_Comm_GsumVecInt(comm, proc_array, itmp, nprocs);*/
   ML_gsum_vec_int(&proc_array, &itmp, nprocs, comm);
   for ( i = nprocs-1; i >= 1; i-- ) proc_array[i] = proc_array[i-1];
   proc_array[0] = 0;
   for ( i = 1; i < nprocs; i++ ) proc_array[i] += proc_array[i-1];
   ML_memory_free( (void**) &itmp );
   (*offset) = proc_array[mypid];
   ML_memory_alloc( (void**) &dtmp, N_total * sizeof(double), "KLC" );
   for ( i = 0; i < N_internal; i++ )
      dtmp[i] = (double) (proc_array[mypid] + i);
   ML_exchange_bdry(dtmp, inmat->comminfo, N_internal, comm, ML_OVERWRITE,NULL);
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

#ifdef out
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
   func->func_ptr( vsolver, 0, NULL, 0, NULL);

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
   func->func_ptr( vsolver, 0, NULL, 0, NULL);
   Amat = (SuperMatrix*) solver->Mat1;
   if (Amat != NULL) {
      Destroy_CompCol_Matrix(Amat);
      ML_memory_free((void**) &Amat);
   }
   solver->Mat1 = NULL;

#else
   solver = (ML_Solver *) vsolver;
   solver->reuse_flag = -999;
   func->func_ptr( vsolver, 0, NULL, 0, NULL);
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
   if ( op->getrow->func_ptr != NULL ) {
      getrow_flag = 1;
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
	flag = op->getrow->func_ptr((void*)op, 1, &i, space-nz_ptr,
                              &(cols[nz_ptr]), &(vals[nz_ptr]), &length);

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
      ML_Set_Amatrix_Matvec(subml, 0, CSR_matvec);
      ML_CommInfoOP_Set_neighbors(&(subml->Amat[0].getrow->pre_comm), 0,
                               NULL, ML_OVERWRITE, NULL, 0);
      ML_Operator_Set_Getrow(&(subml->Amat[0]), 
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
   if ( sl->csolve->func->func_ptr == ML_SuperLU_Solve ) reuse = 1;
   else
   {
      reuse = 0;
      sl->csolve->func->func_ptr = ML_SuperLU_Solve;
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
   if ( op->getrow->func_ptr != NULL ) {
      getrow_flag = 1;
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
	flag = op->getrow->func_ptr((void*)op, 1, &i, space-nz_ptr,
                              &(cols[nz_ptr]), &(vals[nz_ptr]), &length);
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
      ML_Operator_Set_Getrow(&(subml->Amat[0]), 
                             subml->Amat[0].outvec_leng, CSR_getrow);

      ML_Set_Amatrix_Matvec(subml, 0, CSR_matvec);
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
   if ( sl->csolve->func->func_ptr == ML_SuperLU_Solve ) reuse = 1;
   else
   {
      reuse = 0;
      sl->csolve->func->func_ptr = ML_SuperLU_Solve;
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
#endif
#endif

   return 0;
}
#endif

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
#else
   if (ml != NULL) ML_avoid_unused_param((void *) ml);
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
                          (void *) xsfer_op,
                          xsfer_op->Nlocal_rows,
                          ML_OperatorAGX_Restrict, 1);

   ML_Operator_Set_Getrow(&(ml->Rmat[L1]), 
	    (xsfer_op->Nlocal_rows + xsfer_op->Nremote_rows) *stride,
			  ML_OperatorAGX_Getrows);
   ml->Rmat[L1].data_destroy = ML_Operator2AGX_Destroy;

   ML_Operator_Set_1Levels(&(ml->Pmat[L2]), &(ml->SingleLevel[L2]),
                           &(ml->SingleLevel[L1]));
   ML_Operator_Set_ApplyFuncData(&(ml->Pmat[L2]), leng2, leng*stride,
                      (void *) xsfer_op, leng,
                      ML_OperatorAGX_Prolongate, 0);
   ML_Operator_Set_Getrow(&(ml->Pmat[L2]), 
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

/*****************************************************************************/
/* function to partition local subdomain into blocks using metis             */
/*-------------------------------------------------------------------------- */

int ML_Gen_Blocks_Metis(ML *ml, int level, int *nblocks, int **block_list)
{
   *block_list = (int *) ML_allocate(ml->Amat[level].outvec_leng*sizeof(int));
   if (*block_list == NULL)
      pr_error("ML_Gen_Blocks_Metis: out of space\n");
 
   ML_Operator_BlockPartition(&(ml->Amat[level]), ml->Amat[level].outvec_leng,
                             nblocks, *block_list, NULL, NULL, 0,
			      ML_USEMETIS);
   return 0;
}

/*****************************************************************************/
/* Generate a coarse grid matrix suitable for solution with aggregation      */
/* ------------------------------------------------------------------------- */

int ML_Gen_CoarseSolverAggregation(ML *ml_handle, int level, ML_Aggregate *ag)
{
   int            i, j, k, offset, N_local;
   int            reuse, coarsest_level, flag, space, *cols, nz_ptr;
   int            getrow_flag, osize, *row_ptr, length, zero_flag;
   int            local_nlevels, local_clevel;
   double         *vals, dsize, di, *diagonal;
   void           *data;
   ML_1Level      *sl;
   ML_Operator    *op;
   ML_Matrix_DCSR *csr_mat, *csr2_mat;
   ML             *local_ml;
   ML_Aggregate   *newag;
   ML_Solver      *solver;

#ifdef ML_TIMING
   double t0;
   t0 = GetClock();
#endif

   /* ----------------------------------------------------------------- */
   /* extract local matrix using getrow function and store it into a    */
   /* CSR data object                                                   */
   /* ----------------------------------------------------------------- */

   if ( level < 0 || level >= ml_handle->ML_num_levels )
   {
      printf("ML_Gen_CoarseSolverAggregation ERROR : invalid level number.\n");
      exit(-1);
   }
   op      = (ML_Operator *) &ml_handle->Amat[level];
   data    = op->data;
   osize   = op->outvec_leng;
   row_ptr = (int *) ML_allocate(sizeof(int)*(osize+1));
   space   = osize * 5 + 30;
   getrow_flag = 0;
   if      ( op->getrow->func_ptr != NULL ) getrow_flag = 1;
   else
   {
      printf("ML_Gen_CoarseSolverAggregation ERROR : no getrow function.\n");
      exit(-1);
   }

   flag    = 0;

   while (flag == 0)
   {
      cols    = (int    *) ML_allocate(sizeof(int)*space);
      vals    = (double *) ML_allocate(sizeof(double)*space);

      nz_ptr = 0;
      row_ptr[0] = nz_ptr;
      flag = 1;
      for (i = 0; i < osize; i++)
      {
	flag = op->getrow->func_ptr((void*)op, 1, &i, space-nz_ptr,
                              &(cols[nz_ptr]), &(vals[nz_ptr]), &length);

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
      if (flag == 0)
      {
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
   /* form an global matrix                                             */
   /* ----------------------------------------------------------------- */

   ML_memory_alloc((void**) &csr2_mat, sizeof(ML_Matrix_DCSR), "DCR");
   ML_Gen_Amatrix_Global( csr_mat, csr2_mat, ml_handle->comm, &offset);
   csr2_mat->comminfo = NULL;
   ML_free(cols);
   ML_free(vals);
   ML_free(row_ptr);
   ML_free(csr_mat);

   /* ----------------------------------------------------------------- */
   /* set Aggregation as solver                                         */
   /* ----------------------------------------------------------------- */

   coarsest_level = level;
   sl = &(ml_handle->SingleLevel[coarsest_level]);
   if ( sl->csolve->func->func_ptr == ML_CSolve_Aggr ) reuse = 1;
   else
   {
      reuse = 0;
      sl->csolve->func->func_ptr = ML_CSolve_Aggr;
   }

   /* ----------------------------------------------------------------- */
   /* free up previous storage                                          */
   /* ----------------------------------------------------------------- */

   if ( sl->csolve->data != NULL && reuse == 1 )
   {
      solver = (ML_Solver *) sl->csolve->data;
      if ( solver->dble_params1 != NULL )
      {
         ML_memory_free(  (void**) &(solver->dble_params1) );
         solver->dble_params1 = NULL;
      }
#ifdef ML_CPP
      local_ml = (ML *)solver->void_params1;
#else
      local_ml = solver->void_params1;
#endif
      if ( reuse == 1 ) ML_Destroy( &local_ml );
      if (solver->Mat1 != NULL )
      {
         ML_Matrix_DCSR_Destroy((ML_Matrix_DCSR *) (solver->Mat1));
         ML_memory_free(  (void**) &(solver->Mat1) );
         solver->Mat1 = NULL;
      }
   }

   /* ----------------------------------------------------------------- */
   /* create new context                                                */
   /* ----------------------------------------------------------------- */

   ML_Solver_Create( &solver );
   sl->csolve->data = (void *) solver;
   solver->reuse_flag = 0;
   ML_memory_alloc( (void **) &vals, 3 * sizeof(double), "KLI" );
   N_local = osize;
   vals[0]  = (double) N_local;
   vals[1]  = (double) offset;
   vals[2]  = (double) csr2_mat->mat_n;
   solver->dble_params1 = (double *) vals;
   solver->Mat1 = (void *) csr2_mat;
   solver->void_params2 = (void *) ml_handle->comm;

   /* ----------------------------------------------------------------- */
   /* create new context                                                */
   /* ----------------------------------------------------------------- */

   local_nlevels = 10;
   ML_Create( &local_ml, local_nlevels );
   ML_Set_OutputLevel(local_ml, 0);
   ML_Set_ResidualOutputFrequency(local_ml, 0);
   ML_Set_Comm_MyRank(local_ml, 0);
   ML_Set_Comm_Nprocs(local_ml, 1);
   N_local = csr2_mat->mat_n;
   cols    = csr2_mat->mat_ja;
   vals    = csr2_mat->mat_a;
   row_ptr = csr2_mat->mat_ia;
   ML_Init_Amatrix(local_ml,local_nlevels-1,N_local,N_local,(void*) csr2_mat);
   ML_Operator_Set_ApplyFunc(&(local_ml->Amat[local_nlevels-1]),
			     ML_Matrix_DCSR_Matvec);

   local_ml->Amat[local_nlevels-1].data_destroy = ( void (*)(void *)) ML_Matrix_DCSR_Destroy;
   local_ml->Amat[local_nlevels-1].N_nonzeros = csr2_mat->mat_ia[N_local];
   ML_Operator_Set_Getrow(&(local_ml->Amat[local_nlevels-1]),
			  local_ml->Amat[local_nlevels-1].outvec_leng,ML_Matrix_DCSR_Getrow);

   diagonal = (double *) ML_allocate(N_local * sizeof(double));
   for ( i = 0; i < N_local; i++ )
   {
      for ( j = row_ptr[i]; j < row_ptr[i+1]; j++ )
      {
         if ( cols[j] == i ) {diagonal[i] = vals[j]; break;}
      }
   }
   ML_Set_Amatrix_Diag( local_ml, local_nlevels-1, N_local, diagonal);
   ML_free( diagonal );
   ML_Aggregate_Create( &newag );
   if (ml_handle->comm->ML_mypid == 0) ML_Aggregate_Set_OutputLevel(newag,1);
   else                                ML_Aggregate_Set_OutputLevel(newag,0);
   ML_Aggregate_Set_CoarsenScheme_Uncoupled( newag );
   if ( ag != NULL )
      ML_Aggregate_Set_Threshold( newag, ag->curr_threshold );
   if ( ag != NULL )
      ML_Aggregate_Set_DampingFactor( newag, ag->smoothP_damping_factor );
   ML_Aggregate_Set_MaxCoarseSize( newag, 10 );
   ML_Aggregate_Set_PSmootherType( newag, 0 );
   local_clevel = ML_Gen_MGHierarchy_UsingAggregation(local_ml,
                         local_nlevels-1, ML_DECREASING, newag);
   local_clevel = local_nlevels - local_clevel;
   for (k = local_nlevels-1; k > local_clevel; k--)
   {
      ML_Gen_Smoother_SymGaussSeidel(local_ml, k, ML_PRESMOOTHER, 2, 1.);
      ML_Gen_Smoother_SymGaussSeidel(local_ml, k, ML_POSTSMOOTHER, 2, 1.);
   }
   ML_Gen_CoarseSolverSuperLU( local_ml, local_clevel );
   ML_Gen_Solver(local_ml, ML_MGV, local_nlevels-1, local_clevel);
   ML_Aggregate_Destroy( &newag );
   solver->void_params1 = (void *) local_ml;

#ifdef ML_TIMING
   sl->csolve->build_time = GetClock() - t0;
   ml_handle->timing->total_build_time += sl->csolve->build_time;
   if ( ml_handle->comm->ML_mypid == 0 )
      printf("Local Aggregation total setup time = %e\n",
         sl->csolve->build_time);
#endif

   return 0;
}

/* ------------------------------------------------------------------------- */
/* Generate the Hiptmair smoother.                                           */
/* ------------------------------------------------------------------------- */

/* block version */

int ML_Gen_Smoother_BlockHiptmair( ML *ml , int nl, int pre_or_post, int ntimes,
			      ML_Operator **Tmat_array,
			      ML_Operator **Tmat_trans_array,
			      ML_Operator *Tmat_bc,
			      void *edge_smoother, void **edge_args,
			      void *nodal_smoother, void **nodal_args, int type)
     /*
			      int (*edge_smoother )(void), void *edge_args[],
			      int (*nodal_smoother)(void), void *nodal_args[])
     */

{
   ML_Sm_BlockHiptmair_Data *data = NULL;
   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int start_level, end_level, i, status = 1;
   int *BClist=NULL, BClength=0;
   ML_BdryPts *ml_bc;
   char str[80];
#ifdef ML_TIMING
   double         t0;
   t0 = GetClock();
#endif

   if (nl == ML_ALL_LEVELS) {start_level = 0; end_level = ml->ML_num_levels-1;
#ifdef ML_TIMING
      printf("Timing is incorrect when ML_ALL_LEVELS is used with Hiptmair\n");
#endif
}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_Smoother_Hiptmair: cannot set smoother on level %d\n",
	         start_level);
      return 1;
   }

   fun = ML_Smoother_BlockHiptmair;

   if (pre_or_post == ML_PRESMOOTHER)
   {
      for (i = start_level; i <= end_level; i++)
	  {
         /* Get list of Dirichlet bc, if any. */
         ml_bc = ml->SingleLevel[i].BCs;
         if (ML_BdryPts_Check_Dirichlet_Grid(ml_bc))
            ML_BdryPts_Get_Dirichlet_Grid_Info(ml_bc,&BClength,&BClist);
         ML_Smoother_Create_BlockHiptmair_Data(&data);
	     ML_Smoother_Gen_BlockHiptmair_Data(&data, &(ml->Amat[i]),
			          Tmat_array[i], Tmat_trans_array[i], Tmat_bc,
                                  BClength, BClist,
             edge_smoother, edge_args, nodal_smoother, nodal_args );
             data->reduced_smoother = type;
	     ml->pre_smoother[i].data_destroy = ML_Smoother_Destroy_BlockHiptmair_Data;
         sprintf(str,"Hiptmair_pre%d",i);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), 
				      (void *) data, fun, ntimes, 1.0, str);
         ml->pre_smoother[i].pre_or_post = ML_TAG_PRESM;
         BClist = NULL; BClength = 0;
#ifdef ML_TIMING
         ml->pre_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->pre_smoother[i].build_time;
#endif
      }
   }
   else if (pre_or_post == ML_POSTSMOOTHER)
   {
      printf("ML_Gen_Smoother_Hiptmair: ML_POSTSMOOTHER isn't done.\n");
      for (i = start_level; i <= end_level; i++)
	  {
             sprintf(str,"Hiptmair_post%d",i);
             status = ML_Smoother_Set(&(ml->post_smoother[i]),
				      (void *) data, fun, ntimes, 1.0, str);
#ifdef ML_TIMING
         ml->post_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->post_smoother[i].build_time;
#endif

      }

   }
   else if (pre_or_post == ML_BOTH)
   {
      for (i = start_level; i <= end_level; i++)
	  {
         /* Get list of Dirichlet bc, if any. */
         ml_bc = ml->SingleLevel[i].BCs;
         if (ML_BdryPts_Check_Dirichlet_Grid(ml_bc))
            ML_BdryPts_Get_Dirichlet_Grid_Info(ml_bc,&BClength,&BClist);
         ML_Smoother_Create_BlockHiptmair_Data(&data);
	     ML_Smoother_Gen_BlockHiptmair_Data(&data, &(ml->Amat[i]),
			          Tmat_array[i], Tmat_trans_array[i], Tmat_bc,
					   BClength, BClist,
edge_smoother, edge_args, nodal_smoother, nodal_args );
	     data->reduced_smoother = type;
	     ml->post_smoother[i].data_destroy =
			                            ML_Smoother_Destroy_BlockHiptmair_Data;
         sprintf(str,"Hiptmair_pre%d",i);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), 
				      (void *) data, fun, ntimes, 1.0, str);
         ml->pre_smoother[i].pre_or_post = ML_TAG_PRESM;
         sprintf(str,"Hiptmair_post%d",i);
         status = ML_Smoother_Set(&(ml->post_smoother[i]), 
				      (void *) data, fun, ntimes, 1.0, str);
         ml->post_smoother[i].pre_or_post = ML_TAG_POSTSM;
         BClist = NULL; BClength = 0;
#ifdef ML_TIMING
         ml->post_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->post_smoother[i].build_time;
#endif
      }
   }
   else return(pr_error("ML_Gen_Smoother_BlockHiptmair: unknown pre_or_post choice\n"));
   return(status);
}


/******************************************************************************/

int ML_Gen_Smoother_Hiptmair( ML *ml , int nl, int pre_or_post, int ntimes,
			      ML_Operator **Tmat_array,
			      ML_Operator **Tmat_trans_array,
			      ML_Operator *Tmat_bc,
			      void *edge_smoother, void **edge_args,
			      void *nodal_smoother, void **nodal_args, int type)
     /*
			      int (*edge_smoother )(void), void *edge_args[],
			      int (*nodal_smoother)(void), void *nodal_args[])
     */

{
   ML_Sm_Hiptmair_Data *data;
   int (*fun)(ML_Smoother *, int, double *, int, double *);
   int start_level, end_level, i, status = 1;
   int *BClist=NULL, BClength=0;
   ML_BdryPts *ml_bc;
   char str[80];
#ifdef ML_TIMING
   double         t0;
   t0 = GetClock();
#endif

   if (nl == ML_ALL_LEVELS) {start_level = 0; end_level = ml->ML_num_levels-1;
#ifdef ML_TIMING
      printf("Timing is incorrect when ML_ALL_LEVELS is used with Hiptmair\n");
#endif
}
   else { start_level = nl; end_level = nl;}
   if (start_level < 0) {
      printf("ML_Gen_Smoother_Hiptmair: cannot set smoother on level %d\n",
	         start_level);
      return 1;
   }

   fun = ML_Smoother_Hiptmair;

   if (pre_or_post == ML_PRESMOOTHER)
   {
      for (i = start_level; i <= end_level; i++)
	  {
         /* Get list of Dirichlet bc, if any. */
         ml_bc = ml->SingleLevel[i].BCs;
         if (ML_BdryPts_Check_Dirichlet_Grid(ml_bc))
            ML_BdryPts_Get_Dirichlet_Grid_Info(ml_bc,&BClength,&BClist);
         ML_Smoother_Create_Hiptmair_Data(&data);
	     ML_Smoother_Gen_Hiptmair_Data(&data, &(ml->Amat[i]),
			          Tmat_array[i], Tmat_trans_array[i], Tmat_bc,
                      BClength, BClist,
edge_smoother, edge_args, nodal_smoother, nodal_args );
	     data->reduced_smoother = type;
	     ml->pre_smoother[i].data_destroy = ML_Smoother_Destroy_Hiptmair_Data;
         sprintf(str,"Hiptmair_pre%d",i);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), 
				      (void *) data, fun, ntimes, 1.0, str);
         ml->pre_smoother[i].pre_or_post = ML_TAG_PRESM;
         BClist = NULL; BClength = 0;
#ifdef ML_TIMING
         ml->pre_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->pre_smoother[i].build_time;
#endif
      }
   }
   else if (pre_or_post == ML_POSTSMOOTHER)
   {
      for (i = start_level; i <= end_level; i++) {
         /* Get list of Dirichlet bc, if any. */
         ml_bc = ml->SingleLevel[i].BCs;
         if (ML_BdryPts_Check_Dirichlet_Grid(ml_bc))
            ML_BdryPts_Get_Dirichlet_Grid_Info(ml_bc,&BClength,&BClist);
         ML_Smoother_Create_Hiptmair_Data(&data);
	 ML_Smoother_Gen_Hiptmair_Data(&data, &(ml->Amat[i]),
				       Tmat_array[i], Tmat_trans_array[i], Tmat_bc,
				       BClength, BClist,
				       edge_smoother, edge_args, nodal_smoother, nodal_args );
	 data->reduced_smoother = type;
	 ml->post_smoother[i].data_destroy = ML_Smoother_Destroy_Hiptmair_Data;
         sprintf(str,"Hiptmair_post%d",i);
         status = ML_Smoother_Set(&(ml->post_smoother[i]), 
				      (void *) data, fun, ntimes, 1.0, str);
         ml->post_smoother[i].pre_or_post = ML_TAG_POSTSM;
         BClist = NULL; BClength = 0;
#ifdef ML_TIMING
         ml->post_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->post_smoother[i].build_time;
#endif
      }
   }
   else if (pre_or_post == ML_BOTH)
   {
      for (i = start_level; i <= end_level; i++)
	  {
         /* Get list of Dirichlet bc, if any. */
         ml_bc = ml->SingleLevel[i].BCs;
         if (ML_BdryPts_Check_Dirichlet_Grid(ml_bc))
            ML_BdryPts_Get_Dirichlet_Grid_Info(ml_bc,&BClength,&BClist);
         ML_Smoother_Create_Hiptmair_Data(&data);
	     ML_Smoother_Gen_Hiptmair_Data(&data, &(ml->Amat[i]),
			          Tmat_array[i], Tmat_trans_array[i], Tmat_bc,
					   BClength, BClist,
edge_smoother, edge_args, nodal_smoother, nodal_args );
	     data->reduced_smoother = type;
	     ml->post_smoother[i].data_destroy =
			                            ML_Smoother_Destroy_Hiptmair_Data;
         sprintf(str,"Hiptmair_pre%d",i);
         status = ML_Smoother_Set(&(ml->pre_smoother[i]), 
				      (void *) data, fun, ntimes, 1.0, str);
         ml->pre_smoother[i].pre_or_post = ML_TAG_PRESM;
         sprintf(str,"Hiptmair_post%d",i);
         status = ML_Smoother_Set(&(ml->post_smoother[i]), 
				      (void *) data, fun, ntimes, 1.0, str);
         ml->post_smoother[i].pre_or_post = ML_TAG_POSTSM;
         BClist = NULL; BClength = 0;
#ifdef ML_TIMING
         ml->post_smoother[i].build_time = GetClock() - t0;
         ml->timing->total_build_time   += ml->post_smoother[i].build_time;
#endif
      }
   }
   else return(pr_error("ML_Gen_Smoother_Hiptmair: unknown pre_or_post choice\n"));
   return(status);
}

#include "ml_amesos.h"
#define newrap       /* Should always be defined for better performance */

int ML_build_ggb(ML *ml, void *data)
{
  
  ML                    *ml_ggb=NULL;
  int                    Nrows, Ncols, Nlocal, Nnz;
  ML_Operator           *Pmat=NULL, *Qtilde=NULL;
  struct ML_CSR_MSRdata *csr_data, *mydata, *Qtilde_data;
  int                   *NeighborList, *IndList;
  int                    Nneighbors, nprocs, i;
  char str[80];
  double *zdata, *rap, *values;
  int count, j, one=1;
  double *temp;
#ifdef ML_TIMING
   double t0;
   t0 = GetClock();
#endif

  
  mydata   = (struct ML_CSR_MSRdata *) data;
  csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  

  Ncols = mydata->Ncols;
  Nrows = mydata->Nrows;
  Nnz   = mydata->Nnz;
  
  /* Imported information of the Prolongator */
  csr_data->rowptr  =  mydata->rowptr;
  csr_data->columns =  mydata->columns;
  csr_data->values  =  mydata->values;

  if (ml->void_options != NULL) exit(1);  
  ml->void_options = NULL;
  ML_Create( &ml_ggb, 2);

  Pmat = &(ml_ggb->Pmat[0]);

  ML_Operator_halfClone_Init( &(ml_ggb->Amat[1]),
			      &(ml->Amat[ml->ML_finest_level]));

 
  /* Put sizes and function pointers into ml_ggb */
  
  ML_Operator_Set_1Levels(Pmat, &(ml_ggb->SingleLevel[0]),
			  &(ml_ggb->SingleLevel[1]));

  
  /* Only processor 0 has the coarse vector information */

  if (ml_ggb->comm->ML_mypid == 0) 
    ML_Operator_Set_ApplyFuncData(Pmat, Ncols, Nrows, csr_data, Nrows,NULL,0);
   
  else 
    ML_Operator_Set_ApplyFuncData(Pmat, 0, Nrows, csr_data, Nrows, NULL, 0);
  
  
  ML_Operator_Set_Getrow(Pmat, Nrows, CSR_getrow);
  ML_Operator_Set_ApplyFunc (Pmat, CSR_densematvec); 
  /*  ML_Operator_Set_ApplyFunc (Pmat, CSR_matvec);  */

  nprocs  = ml_ggb->comm->ML_nprocs;           /* Number of processors */  
  Nlocal  =  Pmat->invec_leng;                 /* size of coarse grid vector  */
    
  NeighborList =  (int *) ML_allocate(sizeof(int)*(nprocs-1));  
  IndList      =  (int *) ML_allocate(sizeof(int)*Ncols);  
  
  /* processor 0  communicates with all the processors, however all the other
     nodes communicate only with processor 0    */
  
  if (ml_ggb->comm->ML_mypid == 0) {
    
    /* Construct neighbor list for processor 0  */
    for (i = 1; i < nprocs; i++)  NeighborList[i-1] = i;    
    Nneighbors = nprocs-1;
    
    /* indices of in-comming vector (coarse vector)   */
    for (i = 0; i < Ncols; i++) IndList[i] = i;
    
    /* set communication for processor 0                 */
    ML_CommInfoOP_Set_neighbors(&(Pmat->getrow->pre_comm), Nneighbors, 
				NeighborList,ML_OVERWRITE, NULL, 0);
    
    /* information is sent to the other but not recieved  */
    for (i = 0; i < Nneighbors; i++) 
      ML_CommInfoOP_Set_exch_info(Pmat->getrow->pre_comm, NeighborList[i],
				  0,     NULL , Ncols, IndList ); 
    
  }
  else {
    Nneighbors      = 1;
    NeighborList[0] = 0;   
    
    
    /* indices of in-comming vector (coarse vector)  */
    for (i = 0; i < Ncols; i++) IndList[i] = i;
    
    /* set communication for all other processors                 */
    ML_CommInfoOP_Set_neighbors(&(Pmat->getrow->pre_comm), Nneighbors, 
				NeighborList,ML_OVERWRITE, NULL, 0);
    
    /* recieve information from processor 0 */
    ML_CommInfoOP_Set_exch_info(Pmat->getrow->pre_comm, 0,
  		 		Ncols,    IndList , 0, NULL );
 
  }
  Pmat->data_destroy   = ML_OnlyFreeTopLevelDataPtr;


  /*  Pmat->getrow->pre_comm = ML_CommInfoOP_Create();
      sprintf(str,"%d:",ml_ggb->comm->ML_mypid);
      ML_CommInfoOP_Print(Pmat->getrow->pre_comm, str);
  */

  ML_free(NeighborList);
  ML_free(IndList);

    
  /* ML_Operator_Print(Pmat, "Pmat"); */

  ML_Gen_Restrictor_TransP(ml_ggb, 1, 0);
  /* ML_Operator_Set_ApplyFunc (&(ml_ggb->Rmat[1]), CSR_densematvec); */
 


  /* Dense RAP matrix-matrix product */

#ifdef newrap

  zdata = (double *) ML_allocate( Ncols * Nrows*sizeof(double));
  temp  = (double *) ML_allocate( Nrows*sizeof(double));
  values = csr_data->values;

  if (store_AQ == 1) {
    /* Define Qtilde = K*Q;  to be used later in the GGB cycle (GGB first) */
    Qtilde_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
    Qtilde_data->values  =  (double *) ML_allocate((Ncols*Nrows+1)*sizeof(double));
    Qtilde_data->columns =  (int *) ML_allocate(sizeof(int)*(Nrows*Ncols+1));
    Qtilde_data->rowptr  =  (int *) ML_allocate(sizeof(int)*(Nrows+1));
  }

  values = csr_data->values;
  
  /* Compute the A*Q and store in CSR format for later use */
  for (i = 0; i < Ncols; i++) {
    for (j = 0; j < Nrows; j++) temp[j] = values[j*Ncols+i];
    
    ML_Operator_Apply( &(ml_ggb->Amat[1]),Nrows, temp,  Nrows, &(zdata[i*Nrows]));
    
    
    if (store_AQ == 1)
      for (j = 0; j < Nrows; j++) Qtilde_data->values[j*Ncols+i] = zdata[i*Nrows+j];    

}


  if (store_AQ == 1) {
  
    for (j = 0; j < Nrows+1; j++)      Qtilde_data->rowptr[j]  = mydata->rowptr[j];
    for (j = 0; j < Nrows*Ncols; j++)  Qtilde_data->columns[j] = mydata->columns[j];
    
    
    /* Define Qtilde */
    Qtilde = ML_Operator_Create(Pmat->comm);
    ML_CommInfoOP_Clone(  &(Qtilde->getrow->pre_comm),  Pmat->getrow->pre_comm);
    
    
    if (ml_ggb->comm->ML_mypid == 0) 
      ML_Operator_Set_ApplyFuncData(Qtilde, Ncols, Nrows, Qtilde_data,
				    Nrows, NULL, 0);
    
    else 
      
      ML_Operator_Set_ApplyFuncData(Qtilde, 0, Nrows, Qtilde_data,
				    Nrows, NULL, 0);
    
    
    
    ML_Operator_Set_Getrow(Qtilde, Nrows, CSR_getrow);
    ML_Operator_Set_ApplyFunc (Qtilde, CSR_densematvec); 
    Qtilde->data_destroy = ML_CSR_MSRdata_Destroy;
    
    
    ml_ggb->void_options = (void *) Qtilde;
    
  }


  /*
  printf("before the dump\n"); fflush(stdout);  
  ML_Operator_Print(Qtilde, "Qtilde_print");
  ML_CommInfoOP_Print( (Qtilde->getrow->pre_comm),  "Qcomm");
  ML_Operator_Dump(Qtilde, NULL, NULL,"Qtilde", ML_TRUE);
  printf("after the dump\n"); fflush(stdout);  
  while(1 == 1) ; 
  */

  csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  rap = (double *) ML_allocate( Ncols * Ncols * sizeof(double));
  csr_data->columns =  (int *) ML_allocate(sizeof(int)*(Ncols*Ncols+1));
  csr_data->rowptr  =  (int *) ML_allocate(sizeof(int)*(Ncols+1));
 

  count = 0;
  for (i = 0; i < Ncols; i++) {
    for (j = 0; j < Nrows; j++) temp[j] = values[j*Ncols+i];
    
    for (j = 0; j < Ncols; j++) {
      csr_data->columns[count] = j;
      csr_data->rowptr[i] = Ncols*i;

      /*
	rap[count++] = ML_gdot(Nrows, temp,
	&(zdata[j*Nrows]), ml->comm);
      */
	rap[count++] =  DDOT_F77(&Nrows, temp, &one, &(zdata[j*Nrows]), &one);
      
    }    
        
  }

  /* Communicate the local dot products to get the sum */ 
  double *workspace;
  workspace = (double *) ML_allocate( Ncols*Ncols * sizeof(double));
  ML_gsum_vec_double(&rap, &workspace, Ncols*Ncols, ml_ggb->comm); 
  ML_free(workspace);
  
  
  /* Put the RAP matrix result in CSR format */
  csr_data->rowptr[Ncols] = Ncols*Ncols;
  csr_data->values  =  rap;

  if (ml_ggb->comm->ML_mypid == 0)
    ML_Operator_Set_ApplyFuncData( &(ml_ggb->Amat[0]), Ncols, Ncols, csr_data,
				   Ncols, NULL, 0);
  else 
    ML_Operator_Set_ApplyFuncData( &(ml_ggb->Amat[0]), 0, 0, csr_data,
				   0, NULL, 0);
  
  
  ML_Operator_Set_Getrow(&(ml_ggb->Amat[0]), Ncols, CSR_getrow);
  ML_Operator_Set_ApplyFunc (&(ml_ggb->Amat[0]), CSR_densematvec); 
  ml_ggb->Amat[0].getrow->pre_comm = ML_CommInfoOP_Create();
  ml_ggb->Amat[0].data_destroy = ML_CSR_MSRdata_Destroy;


  /* debug the rap */
  /*
    printf("before the dump\n"); fflush(stdout); 
    ML_Operator_Print(&(ml_ggb->Amat[0]), "RAP_print");
  */
    /*
      ML_CommInfoOP_Print( ml_ggb->Amat[0].getrow->pre_comm,  "RAPcomm");
      ML_Operator_Dump(&(ml_ggb->Amat[0]), NULL, NULL,"RAP", ML_TRUE);
      printf("after the dump\n"); fflush(stdout);      
    */
    
  
  /* Free work arrays */
  ML_free(zdata);
  ML_free(temp);


#else
  ML_Gen_AmatrixRAP(ml_ggb, 1, 0);
#endif



#ifdef ML_TIMING
   ml_ggb->Amat[0].build_time = GetClock() - t0;
   ml_ggb->timing->total_build_time   +=    ml_ggb->Amat[0].build_time;
#endif


   switch( ML_ggb_CoarseSolver ) {
   case 1:
     /* superlu, default solver */
     ML_Gen_CoarseSolverSuperLU( ml_ggb, 0);
     break;
   case 2:
     /* amesos, now default is KLU */
     ML_Gen_Smoother_Amesos(ml_ggb, 0, ML_AMESOS_KLU, -1);
     break;
   default:
     printf("ERROR: coarse solver for GGB not correct\n");
     exit( EXIT_FAILURE );
  }

    
  ML_Gen_Solver(ml_ggb, ML_MGV, 1, 0);


  ml->void_options = (void *) ml_ggb;
  
  return 1;
}



void ML_build_ggb_cheap(ML *ml, void *data)
{
  ML                    *ml_ggb;
  int                    Nrows, Ncols, Nlocal, Nnz;
  ML_Operator           *Pmat,  *Qtilde, *Amat_trans;
  struct ML_CSR_MSRdata *csr_data, *mydata;
  int                   *NeighborList, *IndList;
  int                    Nneighbors, nprocs, i;
  char str[80];
  double *zdata, *rap, *values;
  int count, j;
  double *temp;
#ifdef ML_TIMING
   double t0;
   t0 = GetClock();
#endif

  
  mydata   = (struct ML_CSR_MSRdata *) data;
  csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  

  Ncols = mydata->Ncols;
  Nrows = mydata->Nrows;
  Nnz   = mydata->Nnz;


  /* Imported information of the Prolongator */
  csr_data->rowptr  =  mydata->rowptr;
  csr_data->columns =  mydata->columns;
  csr_data->values  =  mydata->values;



  ML_Create( &ml_ggb, 2);

  Pmat = &(ml_ggb->Pmat[0]);

  ML_Operator_halfClone_Init( &(ml_ggb->Amat[1]),
			      &(ml->Amat[ml->ML_finest_level]));

 
  /* Put sizes and function pointers into ml_ggb */
  
  ML_Operator_Set_1Levels(Pmat, &(ml_ggb->SingleLevel[0]),
			  &(ml_ggb->SingleLevel[1]));

  
  /* Only processor 0 has the coarse vector information */

  if (ml_ggb->comm->ML_mypid == 0) 
    ML_Operator_Set_ApplyFuncData(Pmat, Ncols, Nrows, csr_data,
				  Nrows, NULL, 0);
   
  else 
    ML_Operator_Set_ApplyFuncData(Pmat, 0, Nrows, csr_data,
				  Nrows, NULL, 0);
  
  
  ML_Operator_Set_Getrow(Pmat, Nrows, CSR_getrow);
  ML_Operator_Set_ApplyFunc (Pmat, CSR_densematvec); 
  /*  ML_Operator_Set_ApplyFunc (Pmat, CSR_matvec);  */

  nprocs  = ml_ggb->comm->ML_nprocs;           /* Number of processors */  
  Nlocal  =  Pmat->invec_leng;                 /* size of coarse grid vector  */
    
  NeighborList =  (int *) ML_allocate(sizeof(int)*(nprocs-1));  
  IndList      =  (int *) ML_allocate(sizeof(int)*Ncols);  
  
  /* processor 0  communicates with all the processors, however all the other
     nodes communicate only with processor 0    */
  
  if (ml_ggb->comm->ML_mypid == 0) {
    
    /* Construct neighbor list for processor 0  */
    for (i = 1; i < nprocs; i++)  NeighborList[i-1] = i;    
    Nneighbors = nprocs-1;
    
    /* indices of in-comming vector (coarse vector)   */
    for (i = 0; i < Ncols; i++) IndList[i] = i;
    
    /* set communication for processor 0                 */
    ML_CommInfoOP_Set_neighbors(&(Pmat->getrow->pre_comm), Nneighbors, 
				NeighborList,ML_OVERWRITE, NULL, 0);
    
    /* information is sent to the other but not recieved  */
    for (i = 0; i < Nneighbors; i++) 
      ML_CommInfoOP_Set_exch_info(Pmat->getrow->pre_comm, NeighborList[i],
				  0,     NULL , Ncols, IndList ); 
    
  }
  else {
    Nneighbors      = 1;
    NeighborList[0] = 0;   
    
    
    /* indices of in-comming vector (coarse vector)  */
    for (i = 0; i < Ncols; i++) IndList[i] = i;
    
    /* set communication for all other processors                 */
    ML_CommInfoOP_Set_neighbors(&(Pmat->getrow->pre_comm), Nneighbors, 
				NeighborList,ML_OVERWRITE, NULL, 0);
    
    /* recieve information from processor 0 */
    ML_CommInfoOP_Set_exch_info(Pmat->getrow->pre_comm, 0,
  		 		Ncols,    IndList , 0, NULL );
 
  }
  
  /*  Pmat->getrow->pre_comm = ML_CommInfoOP_Create();
      sprintf(str,"%d:",ml_ggb->comm->ML_mypid);
      ML_CommInfoOP_Print(Pmat->getrow->pre_comm, str);
  */
  
  ML_free(NeighborList);
  ML_free(IndList);
  
    
  /* ML_Operator_Print(Pmat, "Pmat"); */

  ML_Gen_Restrictor_TransP(ml_ggb, 1, 0);

  /* Define Qtilde */

  Qtilde = ML_Operator_Create(Pmat->comm);

  ML_2matmult(&(ml_ggb->Rmat[1]), &(ml_ggb->Amat[1]),
	      Qtilde, ML_CSR_MATRIX );
  
  ML_2matmult(Qtilde, &(ml_ggb->Pmat[0]), &(ml_ggb->Amat[0]),
	      ML_MSR_MATRIX );
  
  ML_Operator_Set_ApplyFunc (Qtilde, CSR_densematvec);
  
  
  


#ifdef ML_TIMING
  ml_ggb->Amat[0].build_time = GetClock() - t0;
  ml_ggb->timing->total_build_time   +=    ml_ggb->Amat[0].build_time;
#endif

  ML_Gen_CoarseSolverSuperLU( ml_ggb, 0);
  ML_Gen_Solver(ml_ggb, ML_MGV, 1, 0);
  
  ml_ggb->void_options = (void *) Qtilde;
  ml->void_options     = (void *) ml_ggb;

 }




void ML_build_ggb_fat(ML *ml, void *data)
{
  ML *ml_ggb;
  ML_Aggregate *ag;
  int  Nrows, Ncols, N_levels, num_PDE_eqns ;

  
  ML_Operator *Amat;

  struct ML_CSR_MSRdata *mydata;
  mydata   = (struct ML_CSR_MSRdata *) data;
  
  Amat = &(ml->Amat[ml->ML_finest_level]);
  
  /* coarsest_level = 1; */
  N_levels       = 2;
  num_PDE_eqns   = Amat->num_PDEs;  
  Ncols = mydata->Ncols;
  Nrows = mydata->Nrows;

  /*  for (i = 0; i < Nrows; i++) printf("eig[%d]=%f \n",i,mydata->values[i]); */
  /*  for (i = 0; i < Nrows; i++) printf("Nrows = %d, Ncols = %d \n",Nrows,Ncols); */

  
  printf("\n ========>       Fattening of Prolongation       <========\n\n");
  

  ML_Create( &ml_ggb, N_levels );  /* Haim: make the '2' N_levels */

  ML_Operator_halfClone_Init( &(ml_ggb->Amat[N_levels-1]), /* make 1 N_levels-1 */
			      &(ml->Amat[ml->ML_finest_level]));
			      
  ML_Aggregate_Create( &ag );
  ML_Aggregate_Set_DampingFactor(ag,0.0);
  ML_Aggregate_Set_Threshold(ag, 0.0);
  ML_Aggregate_Set_MaxCoarseSize( ag, 300);
  ML_Aggregate_Set_NullSpace(ag, num_PDE_eqns, Ncols, mydata->values, Nrows);
  ML_Gen_MGHierarchy_UsingAggregation(ml_ggb, N_levels-1,ML_DECREASING, ag);

  /*  ML_Operator_Print(&(ml_ggb->Pmat[0]), "Pmat");  */

  ML_Gen_CoarseSolverSuperLU( ml_ggb, 0);
  ML_Gen_Solver(ml_ggb, ML_MGV, 1, 0);
   
  ml->void_options = (void *) ml_ggb;

 
}







#ifdef WKC

/* WKC -- herein lies all modified code for using Epetra_MultiVectors! */

#include <iostream>
#include <stdio.h>
#include <unistd.h>

int ML_Solve_MGV( ML *ml , const Epetra_MultiVector &in, Epetra_MultiVector &out)
{
/* Copy input to maintain const */
   Epetra_MultiVector  in_temp ( in );
   int    i, leng, dir_leng, *dir_list, k, level;
   double *diag, *scales;


   /* ------------------------------------------------------------ */
   /* initially set the solution to be all 0           	           */
   /* ------------------------------------------------------------ */

   level = ml->ML_finest_level;
   leng = ml->Amat[level].outvec_leng;
   out.PutScalar(0.0);


   /* ------------------------------------------------------------ */
   /* on the fixed boundaries, set them to be each to the solution */
   /* ------------------------------------------------------------ */

   ML_BdryPts_Get_Dirichlet_Eqn_Info(&(ml->BCs[level]),&dir_leng,&dir_list);
   if ( dir_leng != 0 )
   {
      if (ml->Amat[level].diagonal != NULL) {
         ML_DVector_GetDataPtr(ml->Amat[level].diagonal,&diag);
         for ( i = 0; i < dir_leng; i++ ) {
            for ( int KK = 0 ; KK != in.NumVectors() ; KK++ ) {
               k = dir_list[i];
               out[KK][k] = in_temp[KK][k] / diag[k];
            }
         }
      } else {
         diag = NULL;
         for ( i = 0; i < dir_leng; i++ ) {
            for ( int KK = 0 ; KK != in.NumVectors() ; KK++ ) {
               k = dir_list[i];
               out[KK][k] = in_temp[KK][k];
            }
         }
      }
   }

   /* ------------------------------------------------------------ */
   /* normalization                                                */
   /* ------------------------------------------------------------ */

   ML_DVector_GetDataPtr(&(ml->Amat_Normalization[level]), &scales) ;

/* watch out for this !!!!! */
scales = NULL;

/* This can be made much faster by not copying the data above */
   if ( scales != NULL ) {
      for ( i = 0; i < leng; i++ )
      for ( int KK = 0 ; KK != in.NumVectors() ; KK++ )
         in_temp[KK][i] = in[KK][i] / scales[i];
   } else {
      scales = NULL;
   }


   /* ------------------------------------------------------------ */
   /* call MG v-cycle                                              */
   /* ------------------------------------------------------------ */

   ML_Cycle_MG(&(ml->SingleLevel[ml->ML_finest_level]), out,
                in_temp, ML_ZERO, ml->comm, ML_NO_RES_NORM, ml);

   return 0;
}

/*****************************************************************************/
/* solve using V-cycle multigrid                                             */
/*-------------------------------------------------------------------------- */

/* FOR SIMPLICITY,  #ifdef ML_ANALYSIS ... #endif is removed */

double ML_Cycle_MG(ML_1Level *curr, Epetra_MultiVector &ep_sol,
        Epetra_MultiVector &ep_rhs,
	int approx_all_zeros, ML_Comm *comm, int res_norm_or_not, ML *ml)
{

   int         i, lengc, lengf;
   ML_Operator *Amat, *Rmat;
   ML_Smoother *pre,  *post;
   ML_CSolve   *csolve;
#ifdef RAP_CHECK
   double    norm1, norm2;
#endif


   double      res_norm , *normalscales;

   Amat     = curr->Amat;
   Rmat     = curr->Rmat;
   pre      = curr->pre_smoother;
   post     = curr->post_smoother;
   csolve   = curr->csolve;
   lengf    = Amat->outvec_leng;

   /* ------------------------------------------------------------ */
   /* first do the normalization                                   */
   /* ------------------------------------------------------------ */

   ML_DVector_GetDataPtr(curr->Amat_Normalization, &normalscales) ;
   Epetra_MultiVector ep_rhss ( ep_rhs );

   /* ------------------------------------------------------------ */
   /* smoothing or coarse solve                                    */
   /* ------------------------------------------------------------ */
   if (Rmat->to == NULL) {    /* coarsest grid */

/*if ( comm->ML_mypid == 0 ) cout << ep_sol.MyLength() << " " << ep_sol.GlobalLength() << " " << endl; */

      if ( ML_CSolve_Check( csolve ) == 1 ) {
         ML_CSolve_Apply(csolve, lengf, ep_sol, lengf, ep_rhss);
      } else {
         ML_Smoother_Apply(pre, lengf, ep_sol, lengf, ep_rhss, approx_all_zeros);
         ML_Smoother_Apply(post, lengf, ep_sol, lengf, ep_rhss, ML_NONZERO);
      }
      if (res_norm_or_not == ML_COMPUTE_RES_NORM) {
/* Need to do this better! */
            res_norm = 0.0;
      }
   }
   else {

      Epetra_MultiVector ep_res ( ep_rhss );


      /* --------------------------------------------------------- */
      /* pre-smoothing and compute residual                        */
      /* --------------------------------------------------------- */
      ML_Smoother_Apply(pre, lengf, ep_sol, lengf, ep_rhss, approx_all_zeros);

      if ( ( approx_all_zeros != ML_ZERO ) ||
           ( pre->smoother->func_ptr != NULL) )
      {
         ML_Operator_Apply(Amat, lengf, ep_sol, lengf, ep_res);
         for ( i = 0; i < lengf; i++ )
            for ( int LL = 0 ; LL != ep_res.NumVectors() ; LL++ )
               ep_res[LL][i] = ep_rhss[LL][i] - ep_res[LL][i];
      }
      else {
         for ( i = 0; i < lengf; i++ )
            for ( int LL = 0 ; LL != ep_res.NumVectors() ; LL++ )
               ep_res[LL][i] = ep_rhss[LL][i];
      }


      if (res_norm_or_not == ML_COMPUTE_RES_NORM)
/* need to think of something to do here! */
         res_norm = 0.0;
/*         res_norm = sqrt(ML_gdot(lengf, res, res, comm)); */

      lengc = Rmat->outvec_leng;

/* ADDED COMMUNICATION HERE!!! */
      int tot_size_WKC;
      MPI_Allreduce ( (void *)&lengc , (void *)&tot_size_WKC , 1 , MPI_INT ,
                      MPI_SUM , comm->USR_comm );

/*if ( comm->ML_mypid == 0 ) cout << "Recursing down" << endl; */

      Epetra_Map Map ( tot_size_WKC , lengc , 0 , ep_res.Comm() );
      Epetra_MultiVector ep_rhs2 ( Map , ep_res.NumVectors() );
      Epetra_MultiVector ep_sol2 ( Map , ep_res.NumVectors() );

      /* --------------------------------------------------------- */
      /* normalization                                             */
      /* --------------------------------------------------------- */

      ML_DVector_GetDataPtr(curr->Amat_Normalization, &normalscales) ;
      if ( normalscales != NULL )
         for ( i = 0; i < lengf; i++ )
            for ( int KK = 0 ; KK != ep_res.NumVectors() ; KK++ )
               ep_res[KK][i] /= normalscales[i];


      /* ------------------------------------------------------------ */
      /* transform the data from equation to grid space, do grid      */
      /* transfer and then transfer back to equation space            */
      /* ------------------------------------------------------------ */
      if ( ML_Mapper_Check(curr->eqn2grid) == 1 )
      {
         Epetra_MultiVector dTmp ( ep_rhs2 );
         ML_Mapper_Apply(curr->eqn2grid, ep_res, dTmp );
         for ( i = 0; i < lengf; i++ )
            for ( int KK = 0 ; KK != ep_res.NumVectors() ; KK++ )
               ep_res[KK][i] = dTmp[KK][i];
      }

      ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, ep_res, lengc, ep_rhs2);
      if ( ML_Mapper_Check(Rmat->to->grid2eqn) == 1 )
      {
         Epetra_MultiVector dTmp ( ep_rhs2 );
         ML_Mapper_Apply(Rmat->to->grid2eqn, ep_rhs2, dTmp );
         for ( i = 0; i < lengc; i++ )
            for ( int KK = 0 ; KK != ep_sol2.NumVectors() ; KK++ )
                ep_rhs2[KK][i] = dTmp[KK][i];
      }

      ML_DVector_GetDataPtr(Rmat->to->Amat_Normalization,&normalscales);
      if ( normalscales != NULL )
         for ( i = 0; i < lengc; i++ )
            for ( int KK = 0 ; KK != ep_rhs2.NumVectors() ; KK++ )
               ep_rhs2[KK][i] = ep_rhs2[KK][i] * normalscales[i];

      /* --------------------------------------------------------- */
      /* process the next level and transfer back to this level    */
      /* --------------------------------------------------------- */
      ML_Cycle_MG( Rmat->to, ep_sol2, ep_rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if (ml->ML_scheme == ML_MGW)
	ML_Cycle_MG( Rmat->to, ep_sol2, ep_rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);

/*if ( comm->ML_mypid == 0 ) cout << "Recursing up" << endl; */
      /* ------------------------------------------------------------ */
      /* transform the data from equation to grid space, do grid      */
      /* transfer and then transfer back to equation space            */
      /* ------------------------------------------------------------ */
      if ( ML_Mapper_Check(Rmat->to->eqn2grid) == 1 )
      {
         Epetra_MultiVector dTmp ( ep_sol2 );
         ML_Mapper_Apply(Rmat->to->eqn2grid, ep_sol2, dTmp);
         for ( i = 0; i < lengc; i++ )
            for ( int KK = 0 ; KK != ep_sol2.NumVectors() ; KK++ )
                ep_sol2[KK][i] = dTmp[KK][i];
      }

      ML_Operator_ApplyAndResetBdryPts(Rmat->to->Pmat,lengc,ep_sol2,
                                       lengf,ep_res);

      if ( ML_Mapper_Check(curr->grid2eqn) == 1 )
      {
         Epetra_MultiVector dTmp ( ep_res );
         ML_Mapper_Apply(curr->grid2eqn, ep_res, dTmp);
         for ( i = 0; i < lengf; i++ )
            for ( int KK = 0 ; KK != ep_res.NumVectors() ; KK++ )
               ep_res[KK][i] = dTmp[KK][i];
      }

      /* --------------------------------------------------------- */
      /* post-smoothing                                            */
      /* --------------------------------------------------------- */
      for ( i = 0; i < lengf; i++ )
        for ( int KK = 0 ; KK != ep_sol.NumVectors() ; KK++ )
           ep_sol[KK][i] += ep_res[KK][i];

      ML_Smoother_Apply(post, lengf, ep_sol, lengf, ep_rhss, ML_NONZERO);
   }

   return(res_norm);
}

#endif
