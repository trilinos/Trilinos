/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create AMG prolongator                                       */
/* ------------------------------------------------------------------------- */
/* ML_Gen_MGHierarchy_UsingAMG (top level - set up AMG context)              */
/* ML_AMG_Gen_MGHierarchy (2nd level - loop over all levels)                 */
/* ML_AMG_Gen_Prolongator (3rd level - call coarsening routine)              */
/* ML_AMG_Increment_Level (advance to next level - ascending)                */
/* ML_AMG_Decrement_Level (advance to next level - descending)               */
/* ML_AMG_Identity_Getrows (a kludge - to be fixed later)                    */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL)                                       */
/* Date          : October, 2000                                             */
/* ************************************************************************* */
/* ************************************************************************* */

#include <math.h>
#include "ml_amg_genP.h"
#include "ml_op_utils.h"
#include "ml_agg_genP.h"

/* ************************************************************************* */
/* Top level call to RS AMG coarsener - create AMG context and then          */
/* call second level routine to do coarsening based on traversal scheme      */
/* ------------------------------------------------------------------------- */

int ML_Gen_MGHierarchy_UsingAMG(ML *ml, int start, 
                                int increment_or_decrement, ML_AMG *amg)
{
   int    i, j, level, idata=0, nrows, blksize;
   double dnnz;
   ML_AMG *ml_amg;
#ifdef ML_TIMING
   double t0;
#endif

   /* ----------------------------------------------------------------- */
   /* if user does not provide a ML_AMG object, create a default        */
   /* ----------------------------------------------------------------- */

   if ( amg == NULL ) ML_AMG_Create( &ml_amg );
   else               ml_amg = amg;
   ML_AMG_Set_MaxLevels( ml_amg, ml->ML_num_levels);
   ML_AMG_Set_StartLevel( ml_amg, start );

   /* ----------------------------------------------------------------- */
   /* if system AMG is requested, set appropriate parameters            */
   /* ----------------------------------------------------------------- */

   blksize = ml_amg->num_PDE_eqns;
   if ( blksize > 1 && ml_amg->amg_scheme == ML_AMG_SYSTEM_UNKNOWN )
   {
      nrows = ml->Amat[start].outvec_leng;
      if ( nrows % blksize != 0 )
      {
         printf("Gen_AMG ERROR : local nrows not divisible by blksize\n");
         exit(1);
      }
      ML_memory_alloc((void **) &(ml_amg->blk_info), nrows*sizeof(int), "AM1");
      for ( i = 0; i < nrows; i+= blksize )
         for ( j = 0; j < blksize; j++ ) ml_amg->blk_info[i+j] = j;
   }

   /* ----------------------------------------------------------------- */
   /* create multilevel hierarchy                                       */
   /* ----------------------------------------------------------------- */

   idata = ML_gmax_int(idata, ml->comm);
   if ( ml->comm->ML_mypid == 0 && ml_amg->print_flag < ML_Get_PrintLevel())
      ML_AMG_Print(ml_amg);
#ifdef ML_TIMING
   t0 = GetClock();
#endif
   idata = ML_gmax_int(idata, ml->comm);

   if (increment_or_decrement == ML_INCREASING)
   {
      level = ML_AMG_Gen_MGHierarchy(ml, start, ML_AMG_Increment_Level,
                     ML_AMG_Gen_Prolongator, NULL, ml_amg);
   }
   else if (increment_or_decrement == ML_DECREASING)
   {
      level = ML_AMG_Gen_MGHierarchy(ml, start, ML_AMG_Decrement_Level,
                     ML_AMG_Gen_Prolongator, NULL, ml_amg);
   }
   else 
   {
      if ( ml->comm->ML_mypid == 0 ) 
         printf("ML_Gen_MGHierarchy_UsingAMG : unknown inc_or_dec choice\n");
      exit(1);
   }
#ifdef ML_TIMING
   t0 = GetClock() - t0;
   if ( ml->comm->ML_mypid == 0 && ml_amg->print_flag < ML_Get_PrintLevel())
      printf("AMG total setup time = %e\n", t0);
#endif

   /* ----------------------------------------------------------------- */
   /* compute operator complexity                                       */
   /* ----------------------------------------------------------------- */

   if (increment_or_decrement == ML_INCREASING)
      dnnz = (double) ml->Amat[level-start-1].N_nonzeros;
   else if (increment_or_decrement == ML_DECREASING)
      dnnz = (double) ml->Amat[start+1-level].N_nonzeros;
   dnnz = ML_gsum_double( dnnz, ml->comm );
   ml_amg->operator_complexity += dnnz;

   idata = ML_gmax_int(idata, ml->comm);
   if ( ml->comm->ML_mypid == 0 && ml_amg->print_flag < ML_Get_PrintLevel())
      ML_AMG_Print_Complexity(ml_amg);
   idata = ML_gmax_int(idata, ml->comm);

   if ( amg == NULL ) ML_AMG_Destroy( &ml_amg );
   return(level);
}

/* ************************************************************************* */
/* generate multilevel hierarchy given a subroutine for generating           */
/* prolongation operators                                                    */
/* ------------------------------------------------------------------------- */

int ML_AMG_Gen_MGHierarchy(ML *ml, int fine_level,
        int (*next_level)(ML *, int, ML_Operator *, ML_AMG *amg2),
        int (*user_gen_prolongator)(ML *, int, int, void *, ML_AMG *),
        void *data, ML_AMG *amg)
{
   int  level, next, flag, count=1;
#ifdef ML_DEBUG_AMG
   char string[40];
#endif
#ifdef ML_TIMING
   double t0;
#endif

   level = fine_level;
   next  = next_level(ml, level, &(ml->Amat[fine_level]), amg);

   while (next >= 0) 
   {
      flag = user_gen_prolongator(ml,level,next,(void*)&(ml->Amat[level]),amg);
      if (flag < 0) break;

      ML_Gen_Restrictor_TransP(ml, level, next);

      if ( ml->comm->ML_mypid == 0 && amg->print_flag < ML_Get_PrintLevel())
         printf("ML_AMG : generate Galerkin coarse matrix \n");

#ifdef ML_TIMING
      t0 = GetClock();
#endif
      ML_Gen_AmatrixRAP(ml, level, next);
#ifdef ML_TIMING
      t0 = GetClock() - t0;
      if ( ml->comm->ML_mypid == 0 && amg->print_flag < ML_Get_PrintLevel())
         printf("AMG RAP time at level %3d = %e\n", level, t0);
#endif

      if ( ml->comm->ML_mypid == 0 && amg->print_flag < ML_Get_PrintLevel())
      {
         printf("ML_AMG : coarse matrix generated \n");
         printf("-----------------------------------------------\n");
      }
      level = next;
#ifdef ML_DEBUG_AMG
      if ( level == 1 )
      {
         sprintf(string,"AC%d_%d", ml->comm->ML_mypid, level);
         ML_Operator_Print( &(ml->Amat[level]), string);
      }
#endif
      next  = next_level(ml, next, &(ml->Amat[next]), amg);

      count++;
   }
   return(count);
}

/* ************************************************************************* */
/* call coarsening and create prolongation operator                          */
/* ------------------------------------------------------------------------- */

int ML_AMG_Gen_Prolongator(ML *ml,int level, int clevel, void *data,
                           ML_AMG *amg)
{
   int         Ncoarse, Nfine, gNfine, gNcoarse;
   ML_Operator *Amat, *Pmatrix, *AMGIdentity;

#ifdef ML_TIMING
   double t0;
   t0 =  GetClock();
#endif

   Amat     = (ML_Operator *) data;
   Nfine    = Amat->outvec_leng;
   gNfine   = ML_Comm_GsumInt( ml->comm, Nfine);
   ML_AMG_Set_CurrentLevel( amg, level );
   if ( ml->comm->ML_mypid == 0 && amg->print_flag  < ML_Get_PrintLevel())
      printf("+++++++++++++++++++++++++++++++++++++++++++++++\n");
   Pmatrix = ML_Operator_Create(ml->comm);
   Ncoarse  = ML_AMG_Coarsen(amg, Amat, &Pmatrix, ml->comm);
   gNcoarse = ML_Comm_GsumInt( ml->comm, Ncoarse);
   if ( ml->comm->ML_mypid == 0 && amg->print_flag < ML_Get_PrintLevel()) 
      printf("AMG at level %2d = %d\n", level, gNcoarse);
   if ( gNcoarse == 0 || (1.0*gNfine)/(1.0*gNcoarse+0.1) < 1.05 )
   {
      ML_Operator_Destroy(&Pmatrix);
      return -1;
   }
   AMGIdentity = ML_Operator_Create(ml->comm);
   ML_Operator_Set_ApplyFuncData(AMGIdentity, Amat->invec_leng,
                                  Amat->outvec_leng, (void*) Amat,
                                  Amat->matvec->Nrows, NULL, 0);
   ML_Operator_Set_Getrow(AMGIdentity, Amat->getrow->Nrows, 
                          ML_AMG_Identity_Getrows);
   ML_CommInfoOP_Clone(&(AMGIdentity->getrow->pre_comm),Amat->getrow->pre_comm);
   ML_2matmult(AMGIdentity, Pmatrix, &(ml->Pmat[clevel]), ML_CSR_MATRIX );
   ML_Operator_Destroy(&AMGIdentity);
   ML_Operator_Destroy(&Pmatrix);
   ML_Operator_Set_1Levels(&(ml->Pmat[clevel]),
              &(ml->SingleLevel[clevel]), &(ml->SingleLevel[level]));

#ifdef ML_TIMING
   ml->Pmat[clevel].build_time =  GetClock() - t0;
   ml->timing->total_build_time += ml->Pmat[clevel].build_time;
#endif
   return 0;
}

/* ************************************************************************* */
/* function for advancing to the next coarser level with coarse level        */
/* number larger than the fine levels                                        */
/* ------------------------------------------------------------------------- */

int ML_AMG_Increment_Level(ML *ml, int current_level, ML_Operator *Amat,
                           ML_AMG *amg)
{
   int total_size, temp;

   if (current_level == ml->ML_num_levels-1) return(-1);

   total_size = Amat->invec_leng;
   ML_gsum_scalar_int(&total_size, &temp,ml->comm);
   if ( total_size <= amg->max_coarse_size ) return(-1);

   return(current_level+1);
}

/* ************************************************************************* */
/* function for advancing to the next coarser level with coarse level        */
/* number smaller than the fine levels                                       */
/* ------------------------------------------------------------------------- */

int ML_AMG_Decrement_Level(ML *ml, int current_level, ML_Operator *Amat,
                           ML_AMG *amg)
{
   int total_size, temp;

   if (current_level == 0 ) return(-1);

   total_size = Amat->invec_leng;
   ML_gsum_scalar_int(&total_size, &temp, ml->comm);
   if ( total_size <= amg->max_coarse_size ) return(-1);

   return(current_level-1);
}

/* ************************************************************************* */
/* getrow function for identity matrix                                       */
/* ------------------------------------------------------------------------- */

int ML_AMG_Identity_Getrows(ML_Operator *data, int N_requested_rows, 
           int requested_rows[], int allocated_space, int columns[], 
           double values[], int row_lengths[])
{
   if (N_requested_rows > 1) 
   {
      printf("Too bad. This routine only works with 1 row at a time\n");
      exit(1);
   }
   if ( allocated_space == 0 ) return 0;

   columns[0] = requested_rows[0];
   values[0]  = 1.0;
   row_lengths[0] = 1;

   ML_avoid_unused_param(data);
   return(1);
}

