/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create tentative prolongators                                */
/* ************************************************************************* */
/* Author        : Charles Tong (LLNL)                                       */
/* Date          : December, 1999                                            */
/* ************************************************************************* */
/* ************************************************************************* */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_op_utils.h"
#include "ml_utils.h"
#include "ml_agg_Zoltan.h"
#include "ml_agg_user.h"
#include "ml_agg_VBMETIS.h"
#include "ml_viz_stats.h"


/* ************************************************************************* */
/* variables used for parallel debugging  (Ray)                              */
/* ------------------------------------------------------------------------- */

#ifdef ML_AGGR_PARTEST
extern int **global_mapping = NULL, global_nrows, global_ncoarse;
#endif

/* ************************************************************************* */
/* local defines                                                             */
/* ------------------------------------------------------------------------- */

#define ML_AGGR_READY      -11
#define ML_AGGR_NOTSEL     -12
#define ML_AGGR_SELECTED   -13
#define ML_AGGR_SELECTED2  -14
#define ML_AGGR_BDRY       -15
#define ML_AGGR_MINRANK      1
#define ML_AGGR_MAXLINK      2

/* ************************************************************************* */
/* Constructor                                                               */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Create( ML_Aggregate **ag )
{
   ML_memory_alloc( (void **) ag, sizeof(ML_Aggregate), "AG1" );
   (*ag)->ML_id                      = ML_ID_AGGRE;
   (*ag)->print_flag                 = 1;
   (*ag)->ordering                   = 0;
   (*ag)->min_nodes_per_aggregate    = 2;
   (*ag)->max_neigh_already_selected = 0;
   (*ag)->attach_scheme              = ML_AGGR_MAXLINK;
   (*ag)->coarsen_scheme             = ML_AGGR_UNCOUPLED;
   (*ag)->threshold                  = 0.0;
   (*ag)->smoothP_damping_factor     = 4.0/3.0;
   (*ag)->smoothP_damping_sweeps     = NULL;
   (*ag)->smoothP_type               = 0;  /* point type */
   (*ag)->num_PDE_eqns               = 1;
   (*ag)->nullspace_dim              = 1;
   (*ag)->nullspace_vect             = NULL;
   (*ag)->nullspace_corrupted        = ML_EMPTY;
   (*ag)->keep_agg_information       = 0;
   (*ag)->max_levels                 = 0;
   (*ag)->max_coarse_size            = 32;
   (*ag)->begin_level                = 0;
   (*ag)->cur_level                  = 0;
   (*ag)->aggr_info                  = NULL;
   (*ag)->aggr_count                 = NULL;
   (*ag)->drop_tol_for_smoothing     = 0.0;
   (*ag)->fine_complexity            = 0.0;
   (*ag)->nvblocks                   = 0;
   (*ag)->old_RowOmegas              = 0;
   (*ag)->vblock_info                = NULL;
   (*ag)->operator_complexity        = 0.0;
   (*ag)->keep_P_tentative           = ML_NO;
   (*ag)->smooth_existing_P_tentative = ML_NO;
   (*ag)->P_tentative                = NULL;
   (*ag)->use_transpose              = ML_FALSE;
   (*ag)->Restriction_smoothagg_transpose = ML_FALSE;
/*MS*/
   (*ag)->coarsen_scheme_level       = NULL;
   (*ag)->aggr_options               = NULL;
   (*ag)->aggr_viz_and_stats         = NULL;
   (*ag)->nodal_coord                = NULL;
   (*ag)->field_of_values            = NULL;
/*MS*/
   (*ag)->block_scaled_SA            =  0;
   (*ag)->phase3_agg_creation        = .5;
/*mgee*/
   (*ag)->vblock_data                = NULL;
/*MS*/
   (*ag)->minimizing_energy          = 0;
   (*ag)->minimizing_energy_droptol  = 0.0;
   (*ag)->cheap_minimizing_energy    = 0;
   (*ag)->coarsen_rate               = -1;
   (*ag)->semicoarsen_levels         = -1;
   (*ag)->semicoarsen_coordinate     = 'z';
/*cms*/
   (*ag)->rowsum_threshold           = -1.0; /* defaults to off */

#if defined(AZTEC) && defined(ML_AGGR_READINFO)
   ML_Aggregate_AztecRead(*ag);
#endif

   return 0;
}

/* ************************************************************************* */
/* destructor                                                                */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Destroy( ML_Aggregate **ag )
{
   int i;

   if (*ag != NULL)
   {
      if ( (*ag)->ML_id != ML_ID_AGGRE )
      {
         printf("ML_Aggregate_Destroy : wrong object. \n");
         exit(-1);
      }
      if ((*ag)->nullspace_vect != NULL)
      {
         ML_memory_free((void **)&((*ag)->nullspace_vect));
      }
      if ((*ag)->aggr_info != NULL)
      {
         for ( i = 0; i < (*ag)->max_levels; i++ )
            if ((*ag)->aggr_info[i] != NULL)
               ML_memory_free((void **)&((*ag)->aggr_info[i]));
         ML_memory_free((void **)&((*ag)->aggr_info));
      }
      if ((*ag)->aggr_count != NULL)
      {
         ML_memory_free((void **)&((*ag)->aggr_count));
      }
      if ( (*ag)->P_tentative != NULL)
	ML_Operator_ArrayDestroy( (*ag)->P_tentative, (*ag)->max_levels);

/*MS*/
      if( (*ag)->coarsen_scheme_level != NULL ) {
	ML_memory_free( (void **) &((*ag)->coarsen_scheme_level) );
	(*ag)->coarsen_scheme_level = NULL;
      }
      if( (*ag)->aggr_options != NULL ) {
	ML_memory_free( (void **) &((*ag)->aggr_options) );
	(*ag)->aggr_options = NULL;
      }
      /* aggr_viz_and_stats is cleaned by calling the function
	 `ML_Aggregate_Viz_Stats_Clean', in file "Utils/ml_agg_info.c" */

      if ((*ag)->field_of_values != NULL)
	ML_free((*ag)->field_of_values);
      if ((*ag)->nodal_coord != NULL) {
        /* MS start from 1 because we do not
         * free the finest-level coordinates (given by the
         * user). Recall that in nodal_coord the finest level is
         * always in position 0. */
        for (i = 1 ; i < (*ag)->max_levels ; ++i) {
          if ((*ag)->nodal_coord[i] != NULL)
            ML_free((*ag)->nodal_coord[i]);
        }
	ML_free( (*ag)->nodal_coord );
      }
/*MS*/
/*mgee*/
      if ((*ag)->vblock_data != NULL)
      {
         for ( i = 0; i < (*ag)->max_levels; i++ )
            ML_Aggregate_Destroy_Vblocks_CoarsenScheme_VBMETIS((*ag),i);
         ML_free((*ag)->vblock_data);
         (*ag)->vblock_data = NULL;
      }
      if ((*ag)->smoothP_damping_sweeps != NULL)
        ML_free((*ag)->smoothP_damping_sweeps);

      ML_memory_free( (void **) ag );
      (*ag) = NULL;

   }
   return 0;
}



/* Steers how the MIS  and Uncoupled handle phase 3 of aggregation.
 * Values near 0 create few additional aggregates.Large values create
 * many additional aggregates. Convergence can be improve convergence
 * by new aggregates but nonzero fill-in increases on coarse meshes.
 * Default: .5. The basic idea is that the ratio of vertex neighbors
 * that are already aggregated to the total number of vertex neighbors
 * is compared with alpha^factor where alpha is the ratio of vertices
 * aggregated in phase 1 to the total number of vertices. This is used
 * to see if a new aggregate should be created.
 */
int ML_Aggregate_Set_Phase3AggregateCreationAggressiveness(
			   ML_Aggregate *ag, double factor) {

   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_Phase3AggregateCreationAggressiveness: wrong object. \n");
      exit(-1);
   }
   if (factor < 0) {
      printf("ML_Aggregate_Set_Phase3AggregateCreationAggressiveness: negative settings (%e) not allowed.\n",factor);
      exit(-1);
   }
   ag->phase3_agg_creation        = factor;
   return 0;
}

/* ************************************************************************* */
/* set minimum number of nodes per aggregate for phase 2                     */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_MinNodesPerAggregate( ML_Aggregate *ag, int nnodes )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_MinNodesPerAggregate : wrong object. \n");
      exit(-1);
   }
   if ( nnodes <= 1 ) ag->min_nodes_per_aggregate = 1;
   else               ag->min_nodes_per_aggregate = nnodes;
   return 0;
}
/* ************************************************************************* */
/* set Prolongator smoother to use block diagonal scaling                    */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_BlockDiagScaling( ML_Aggregate *ag)
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_BlockDiagScaling : wrong object. \n");
      exit(-1);
   }
   ag->block_scaled_SA = 1;
   return 0;
}
int ML_Aggregate_KeepInfo(ML_Aggregate *ag, int value)
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_BlockDiagScaling : wrong object. \n");
      exit(-1);
   }
   ag->keep_agg_information = value;
   return 0;
}


/* ************************************************************************* */
/* set Prolongator smoother to use point diagonal scaling                    */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_PointDiagScaling( ML_Aggregate *ag)
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_BlockDiagScaling : wrong object. \n");
      exit(-1);
   }
   ag->block_scaled_SA = 0;
   return 0;
}


/* ************************************************************************* */
/* set output level                                                          */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_OutputLevel( ML_Aggregate *ag, double level )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_OutputLevel : wrong object. \n");
      exit(-1);
   }
   ag->print_flag = level;
   return 0;
}

int ML_Aggregate_Set_Reuse(ML_Aggregate *ag)
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_Reuse : wrong object. \n");
      exit(-1);
   }
   ag->keep_P_tentative = ML_YES;
   return 0;
}

/* ************************************************************************* */
/* set random or natural ordering for traversing the matrix graph            */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_NaturalOrdering( ML_Aggregate *ag )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_NaturalOrdering : wrong object. \n");
      exit(-1);
   }
   ag->ordering = 0;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_RandomOrdering( ML_Aggregate *ag )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_RandomOrdering : wrong object. \n");
      exit(-1);
   }
   ag->ordering = 1;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_GraphOrdering( ML_Aggregate *ag )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_GraphOrdering : wrong object. \n");
      exit(-1);
   }
   ag->ordering = 2;
   return 0;
}

/* ************************************************************************* */
/* select scheme to put un-aggregated nodes into existing aggregates         */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_AttachScheme_MaxLink( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_AttachScheme_MaxLink : wrong object. \n");
      exit(-1);
   }
   ag->attach_scheme = ML_AGGR_MAXLINK;
   return 0;
}

/* -------------------------------------------------------------------- */

int ML_Aggregate_Set_AttachScheme_MinRank( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_AttachScheme_MinRank : wrong object. \n");
      exit(-1);
   }
   ag->attach_scheme = ML_AGGR_MINRANK;
   return 0;
}

/* ************************************************************************* */
/* set maximum coarsest grid size                                            */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_MaxCoarseSize( ML_Aggregate *ag, int size  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_MaxCoarseSize : wrong object. \n");
      exit(-1);
   }
   ag->max_coarse_size = size;
   return 0;
}

/* ************************************************************************* */
/* select coarsening scheme                                             */
/* -------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_Uncoupled( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_Uncoupled : wrong object. \n");
      exit(-1);
   }
   ag->coarsen_scheme = ML_AGGR_UNCOUPLED;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_Coupled( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_Coupled : wrong object. \n");
      exit(-1);
   }
   ag->coarsen_scheme = ML_AGGR_COUPLED;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_MIS( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_MIS : wrong object. \n");
      exit(-1);
   }
   ag->coarsen_scheme = ML_AGGR_MIS;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_UncoupledMIS( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_UncoupledMIS : wrong object. \n");
      exit(-1);
   }
   ag->coarsen_scheme = ML_AGGR_HYBRIDUM;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_UncoupledCoupled( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_UncoupledCoupled:wrong object\n");
      exit(-1);
   }
   ag->coarsen_scheme = ML_AGGR_HYBRIDUC;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_DD( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_DD : wrong object. \n");
      exit(-1);
   }
   ag->coarsen_scheme = ML_AGGR_DD;
   return 0;
}

/*MS*/
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_METIS( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_METIS : wrong object. \n");
      exit(-1);
   }
   ag->coarsen_scheme = ML_AGGR_METIS;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_ParMETIS( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_ParMETIS : wrong object. \n");
      exit(-1);
   }
   ag->coarsen_scheme = ML_AGGR_PARMETIS;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_Zoltan( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_Zoltan : wrong object. \n");
      exit(-1);
   }
   ag->coarsen_scheme = ML_AGGR_ZOLTAN;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenScheme_User( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_Zoltan : wrong object. \n");
      exit(-1);
   }
   ag->coarsen_scheme = ML_AGGR_USER;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Get_CoarsenScheme( ML_Aggregate *ag  )
{
   if ( ag->ML_id != ML_ID_AGGRE )
      pr_error("ML_Aggregate_Get_CoarsenScheme: wrong object. \n");
   return ag->coarsen_scheme;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenSchemeLevel( int level, int MaxLevels,
					 ML_Aggregate *ag,
					 int choice )
{

  int i;

  if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CoarsenScheme_METIS : wrong object. \n");
      exit(-1);
   }
   if( ag->coarsen_scheme_level == NULL ) {
     ML_memory_alloc((void**) &(ag->coarsen_scheme_level),
		     sizeof(int)*MaxLevels,"coarsen_scheme_level");
     for( i=0 ; i<MaxLevels ; i++ )
       ag->coarsen_scheme_level[i] = choice;
   }

   if( level<-1 || level>=MaxLevels ) {
     fprintf( stderr,
	      "*ML*ERR* level not valid (%d), MaxLevels=%d\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      level, MaxLevels,
	      __FILE__,
	      __LINE__ );
     return 1;
   }

   if( level == -1 ) {
     for( i=0 ; i<MaxLevels ; i++ )
       ag->coarsen_scheme_level[i] = choice;
   } else {
     ag->coarsen_scheme_level[level] = choice;
   }

   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenSchemeLevel_Coupled( int level, int MaxLevels,
						 ML_Aggregate *ag  )
{
  return( ML_Aggregate_Set_CoarsenSchemeLevel(level, MaxLevels,
					      ag, ML_AGGR_COUPLED) );
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenSchemeLevel_UncoupledMIS( int level, int MaxLevels,
						 ML_Aggregate *ag  )
{
  return( ML_Aggregate_Set_CoarsenSchemeLevel(level, MaxLevels,
					      ag, ML_AGGR_HYBRIDUM) );
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenSchemeLevel_Uncoupled( int level, int MaxLevels,
						   ML_Aggregate *ag  )
{
  return( ML_Aggregate_Set_CoarsenSchemeLevel(level, MaxLevels,
					      ag, ML_AGGR_UNCOUPLED) );
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenSchemeLevel_MIS( int level, int MaxLevels,
					     ML_Aggregate *ag  )
{
  return( ML_Aggregate_Set_CoarsenSchemeLevel(level, MaxLevels,
					      ag, ML_AGGR_MIS) );
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenSchemeLevel_METIS(int level, int MaxLevels,
					      ML_Aggregate *ag)
{
  return( ML_Aggregate_Set_CoarsenSchemeLevel(level, MaxLevels,
					      ag, ML_AGGR_METIS));
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenSchemeLevel_ParMETIS(int level, int MaxLevels,
						 ML_Aggregate *ag)
{
  return( ML_Aggregate_Set_CoarsenSchemeLevel(level, MaxLevels,
					      ag, ML_AGGR_PARMETIS));
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenSchemeLevel_User(int level, int MaxLevels,
                                             ML_Aggregate *ag)
{
  return( ML_Aggregate_Set_CoarsenSchemeLevel(level, MaxLevels,
					      ag, ML_AGGR_USER));
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenSchemeLevel_Zoltan( int level, int MaxLevels,
						  ML_Aggregate *ag  )
{
  return( ML_Aggregate_Set_CoarsenSchemeLevel(level, MaxLevels,
					      ag, ML_AGGR_ZOLTAN) );
}
/*ms*/

/* ************************************************************************* */
/* set/reset aggregation threshold                                           */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_Threshold( ML_Aggregate *ag, double epsilon )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_Threshold : wrong object. \n");
      exit(-1);
   }
   if ( epsilon > 0.0 ) ag->threshold = epsilon;
   else                 ag->threshold = 0.0;

   ag->drop_tol_for_smoothing = ag->threshold;

   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Reset_Threshold( ML_Aggregate *ag )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Reset_Threshold : wrong object. \n");
      exit(-1);
   }
   ag->threshold = 0.0;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_RowSum_Threshold( ML_Aggregate *ag, double epsilon )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_RowSum_Threshold : wrong object. \n");
      exit(-1);
   }
   if ( epsilon > 0.0 ) ag->rowsum_threshold = epsilon;
   else                 ag->rowsum_threshold = -1.0;

   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Reset_RowSum_Threshold( ML_Aggregate *ag )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Reset_RowSum_Threshold : wrong object. \n");
      exit(-1);
   }
   ag->rowsum_threshold = -1.0;
   return 0;
}

/* ************************************************************************* */
/* Set flag controlling whether to use existing tentative prolongator.       */
/* ************************************************************************* */

int ML_Aggregate_Set_Flag_SmoothExistingTentativeP( ML_Aggregate *ag, int flag)
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_DampingFactor : wrong object. \n");
      exit(-1);
   }
   ag->smooth_existing_P_tentative = flag;
   return 0;
}

/* ************************************************************************* */
/* Retrieve flag controlling whether to use existing tentative prolongator.  */
/* ************************************************************************* */

int ML_Aggregate_Get_Flag_SmoothExistingTentativeP( ML_Aggregate *ag)
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_DampingFactor : wrong object. \n");
      exit(-1);
   }
   return ag->smooth_existing_P_tentative;
}


/* ************************************************************************* */
/* set damping factor for the smoothed prolongator                           */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_DampingFactor( ML_Aggregate *ag, double factor )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_DampingFactor : wrong object. \n");
      exit(-1);
   }
   ag->smoothP_damping_factor = factor;
   return 0;
}


/* ************************************************************************* */
/* set number of damping sweeps for the smoothed prolongator                 */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_DampingSweeps( ML_Aggregate *ag, int numSweeps, int level )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_DampingFactor : wrong object. \n");
      exit(-1);
   }
   if (ag->smoothP_damping_sweeps==NULL)
     pr_error("ML_Aggregate_Set_DampingSweeps:  Memory not allocated.  Call ML_Aggregate_Set_MaxLevels() first.\n");
   if (level == ML_ALL_LEVELS) {
     int i;
     for (i=0; i<ag->max_levels; i++) ag->smoothP_damping_sweeps[i] = numSweeps;
   }
   else
     ag->smoothP_damping_sweeps[level] = numSweeps;
   return 0;
}


/* ************************************************************************* */
/* retrieve number of damping sweeps for the smoothed prolongator            */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Get_DampingSweeps( ML_Aggregate *ag, int level)
{
   if ( ag->ML_id != ML_ID_AGGRE )
     pr_error("ML_Aggregate_Set_DampingFactor : wrong object. \n");
   if (level >= ag->max_levels)
     pr_error("ML_Aggregate_Get_DampingSweeps: largest allowable level = %d\n",
              ag->max_levels);
   return ag->smoothP_damping_sweeps[level];
}

/* ************************************************************************* */
/* set prolongator smoother type (diagonal, block diagonal = 0, 1)      */
/* -------------------------------------------------------------------- */

int ML_Aggregate_Set_PSmootherType( ML_Aggregate *ag, int stype )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_PSmootherType : wrong object. \n");
      exit(-1);
   }
   if ( stype == 0 ) ag->smoothP_type = stype;
   else              ag->smoothP_type = 1;
   return 0;
}

/* ************************************************************************* */
/* set max number of levels for the smoothed prolongator                     */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_MaxLevels( ML_Aggregate *ag, int level )
{
   int i;

   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_MaxLevels : wrong object. \n");
      exit(-1);
   }
   if (ag->max_levels != 0) {
     if (ag->max_levels != level)
       pr_error("ML_Aggregate_Set_MaxLevels : max_levels is already set.\n");
     else
       return 0;
   }
   ag->max_levels = level;
#ifdef ML_CPP
   ML_memory_alloc((void**) &(ag->aggr_info), level*sizeof(int*),"AGu");
#else
   ML_memory_alloc((void*) &(ag->aggr_info), level*sizeof(int*),"AGu");
#endif
   for ( i = 0; i < level; i++ ) ag->aggr_info[i] = NULL;
#ifdef ML_CPP
   ML_memory_alloc((void**) &(ag->aggr_count), level*sizeof(int),"AGx");
#else
   ML_memory_alloc((void*) &(ag->aggr_count), level*sizeof(int),"AGx");
#endif
   if (ag->smoothP_damping_sweeps == NULL) {
     ag->smoothP_damping_sweeps = (int*) ML_allocate(level*sizeof(int));
     ML_Aggregate_Set_DampingSweeps(ag,1,ML_ALL_LEVELS);
   }
   else
     pr_error("ML_Aggregate_Set_MaxLevels: array 'smoothP_damping_sweeps' already allocated\n");
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CurrentLevel( ML_Aggregate *ag, int level )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_CurrentLevel : wrong object. \n");
      exit(-1);
   }
   ag->cur_level = level;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_StartLevel( ML_Aggregate *ag, int level )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Set_StartLevel : wrong object. \n");
      exit(-1);
   }
   ag->begin_level = level;
   return 0;
}

/* ************************************************************************* */
/* get the aggregation information for a given level                         */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Get_AggrCount( ML_Aggregate *ag, int level )
{
   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Get_AggrCount : wrong object. \n");
      exit(-1);
   }
   if ( level < -1 || level >= ag->max_levels )
   {
      printf("ML_Aggregate_Get_AggrCount : level number not valid. \n");
      exit(-1);
   }
   return ag->aggr_count[level];
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Get_AggrMap( ML_Aggregate *ag, int level, int **map )
{

   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Get_AggrMap : wrong object. \n");
      exit(-1);
   }
   if ( level < 0 || level >= ag->max_levels )
   {
      printf("ML_Aggregate_Get_AggrMap : level number not valid. \n");
      exit(-1);
   }
   (*map) = ag->aggr_info[level];
   return 0;
}

/* ************************************************************************* */
/* set nullspace                                                             */
/* null_vect can be NULL, in which case the default nullspace                */
/* (ones for each PDE) will be used.  leng is the dimension of A at the      */
/* coarsest level                                                            */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_NullSpace(ML_Aggregate *ag, int num_PDE_eqns,
                               int null_dim, double *null_vect, int leng)
{
   int nbytes, i;
#ifdef ML_AGGR_NSOUTPUT
   FILE *fp;
   char fname[100];
#endif

   if ((null_vect == NULL) && (num_PDE_eqns != null_dim))
   {
      printf("WARNING:  When no nullspace vector is specified, the number\n");
      printf("of PDE equations must equal the nullspace dimension.\n");
   }

   /* first set the 2 integer arguments */

   ag->num_PDE_eqns  = num_PDE_eqns;
   ag->nullspace_dim = null_dim;

   /* if there was a nullspace vector specified before, free it */
   if (ag->nullspace_vect != NULL) {
     /* if this is the fine grid, indicate that we have overwritten */
     /* the user's nullspace.                                       */
     if (ag->nullspace_corrupted == ML_EMPTY)
       ag->nullspace_corrupted = ML_YES;

     ML_memory_free((void **)&(ag->nullspace_vect));
   }

   /* if the user-supplied nullspace vector isn't null, allocate space */
   /* and load it */

#ifdef ML_AGGR_NSOUTPUT
   sprintf( fname, "null.%d", global_comm->ML_mypid);
   fp = fopen(fname, "w");
#endif

   if (null_vect != NULL)
   {
     /* If the fine grid operator has no null space, indicate */
     /* that the fine grid nullspace has not been corrupted.  */

     if (ag->nullspace_corrupted == ML_EMPTY)
       ag->nullspace_corrupted = ML_NO;

      nbytes = leng * null_dim * sizeof(double);

      ML_memory_alloc((void **)&(ag->nullspace_vect), (unsigned int) nbytes, "ns");

      for (i=0; i < leng*null_dim; i++)
      {
         (ag->nullspace_vect)[i] = null_vect[i];
#ifdef ML_AGGR_NSOUTPUT
         fprintf(fp, "null %5d (%5d) = %e\n", i, leng*null_dim, null_vect[i]);
#endif
      }

   }
   else
      ag->nullspace_vect = NULL;

#ifdef ML_AGGR_NSOUTPUT
   fclose(fp);
#endif
   return 0;
}

/* ************************************************************************* */
/* scale the null space using the vector 'scale_vect'. Note: the nullspace   */
/* can be NULL, in which case the default nullspace (ones for each PDE) will */
/* be used.  length is the local dimension of A                              */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Scale_NullSpace(ML_Aggregate *ag, double *scale_vect,
				 int length)
{
   int nbytes, j, k;
   double *null_vect;
   int num_PDE_eqns, null_dim;

   /* first pull out null space information */

   num_PDE_eqns = ag->num_PDE_eqns;
   null_dim     = ag->nullspace_dim;
   null_vect    = ag->nullspace_vect;


   if ((null_vect == NULL) && (num_PDE_eqns != null_dim))
   {
      printf("WARNING:  When no nullspace vector is specified, the number\n");
      printf("of PDE equations must equal the nullspace dimension.\n");
   }


   /* if the user-supplied nullspace vector is null, allocate space */
   /* and load it */

   if (null_vect == NULL) {
     nbytes = length * null_dim * sizeof(double);

     ML_memory_alloc((void **)&(ag->nullspace_vect), (unsigned int) nbytes, "ns");
     null_vect = ag->nullspace_vect;

     for (j = 0; j < length; j++) {
       for (k = 0; k < null_dim; k++) {
	  if (j % num_PDE_eqns == k) null_vect[k*length+j] = 1.0;
	  else                       null_vect[k*length+j] = 0.0;
       }
      }
   }

   if (scale_vect == NULL) {
     printf("ML_Aggregate_Scale_NullSpace: scale vector is null\n");
     return 1;
   }

   for (k = 0; k < null_dim; k++) {
      for (j = 0; j < length; j++) null_vect[k*length+j] /= scale_vect[j];
   }

   return 0;
}

/* ************************************************************************* */
/* ************************************************************************* */
/* Coarsening routine                                                        */
/* ------------------------------------------------------------------------- */
int ML_Aggregate_Coarsen( ML_Aggregate *ag, ML_Operator *Amatrix,
                          ML_Operator **Pmatrix, ML_Comm *comm)
{
   int i=1, ndofs, Ncoarse, coarsen_scheme, status;
   int mypid, nprocs;
   char *label;

#ifdef ML_TIMING
   double t0;
   t0 = GetClock();
#endif

   label = ML_memory_check(NULL);
   if (label != NULL)
     if ( label[0] == 'L')
       if ( ( label[2] == ':') || ( label[3] == ':') )
	 sscanf(&(label[1]),"%d",&i);

   if (i != 1) ML_memory_check("L%d: agg start",i);
   else ML_memory_check("agg start");

   mypid = comm->ML_mypid;

   if ( ag->ML_id != ML_ID_AGGRE )
   {
      printf("ML_Aggregate_Coarsen : wrong object. \n");
      exit(-1);
   }

   if (mypid == 0 && ag->print_flag < ML_Get_PrintLevel())
      printf("ML_Aggregate_Coarsen (level %d) begins\n", ag->cur_level);

/* #### moved this somewhere else ?? */
   Amatrix->num_PDEs = ag->num_PDE_eqns;
   Amatrix->num_rigid = ag->nullspace_dim;

   /* check to see which aggregation algorithm to use */

   ndofs = Amatrix->outvec_leng;
   /*MS*/
   if( ag->coarsen_scheme_level == NULL ) {
     coarsen_scheme = ag->coarsen_scheme;
   } else {
     coarsen_scheme = ag->coarsen_scheme_level[ag->cur_level];
   }
   /*ms*/
   if (coarsen_scheme == ML_AGGR_HYBRIDUM) {
     if ( ndofs < 250 ) ndofs = 0; else ndofs = 1;
   }
   else {
     if ( ndofs < 2 ) ndofs = 0; else ndofs = 1;
   }
   ML_gsum_scalar_int(&ndofs, &i, comm);
   /*
   gmin = ML_gmax_double((double) (-1.0 * Amatrix->outvec_leng) , comm);
   if (comm->ML_mypid == 0)
      printf("Smallest Amatrix->outvec_leng = %d\n", (int) (-1 * gmin));
   */
   nprocs = comm->ML_nprocs;
   if ( ndofs == nprocs )
   {
      if (coarsen_scheme == ML_AGGR_UNCOUPLED)
         coarsen_scheme = ML_AGGR_UNCOUPLED;
      else if (coarsen_scheme == ML_AGGR_COUPLED)
         coarsen_scheme = ML_AGGR_COUPLED;
      else if (coarsen_scheme == ML_AGGR_MIS)
         coarsen_scheme = ML_AGGR_MIS;
      else if (coarsen_scheme == ML_AGGR_HYBRIDUC)
         coarsen_scheme = ML_AGGR_UNCOUPLED;
      else if (coarsen_scheme == ML_AGGR_HYBRIDUM)
         coarsen_scheme = ML_AGGR_UNCOUPLED;
/*MS*/
      else if (coarsen_scheme == ML_AGGR_METIS)
         coarsen_scheme = ML_AGGR_METIS;
      else if (coarsen_scheme == ML_AGGR_PARMETIS)
	coarsen_scheme = ML_AGGR_PARMETIS;
      else if (coarsen_scheme == ML_AGGR_ZOLTAN)
	coarsen_scheme = ML_AGGR_ZOLTAN;
      else if (coarsen_scheme == ML_AGGR_USER)
	coarsen_scheme = ML_AGGR_USER;
/*ms*/
/*mgee*/
      else if (coarsen_scheme == ML_AGGR_VBMETIS)
         coarsen_scheme = ML_AGGR_VBMETIS;
      else
         coarsen_scheme = ML_AGGR_UNCOUPLED;
   }
   else
   {
/* JJH why????????????????????? */
/*#ifdef ML_repartition*/
      if (coarsen_scheme == ML_AGGR_UNCOUPLED)
         coarsen_scheme = ML_AGGR_UNCOUPLED;
      else
/*#endif*/
      if (coarsen_scheme == ML_AGGR_COUPLED)
         coarsen_scheme = ML_AGGR_COUPLED;
      else if (coarsen_scheme == ML_AGGR_MIS)
         coarsen_scheme = ML_AGGR_MIS;
      else if (coarsen_scheme == ML_AGGR_HYBRIDUC)
         coarsen_scheme = ML_AGGR_COUPLED;
      else if (coarsen_scheme == ML_AGGR_HYBRIDUM)
         coarsen_scheme = ML_AGGR_MIS;
/*MS*/
      else if (coarsen_scheme == ML_AGGR_METIS)
         coarsen_scheme = ML_AGGR_METIS;
      else if (coarsen_scheme == ML_AGGR_PARMETIS)
	coarsen_scheme = ML_AGGR_PARMETIS;
      else if (coarsen_scheme == ML_AGGR_ZOLTAN)
	coarsen_scheme = ML_AGGR_ZOLTAN;
      else if (coarsen_scheme == ML_AGGR_USER)
	coarsen_scheme = ML_AGGR_USER;
/*ms*/
      else
      {
         /* The following can cause a memory leak.  This is now
            freed in ML_AGG_Gen_Prolongator() in ml_agg_genP.c. */
         /*(*Pmatrix) = NULL;*/
         return 0;
      }
   }

   if (coarsen_scheme == ML_AGGR_MIS) {
      i = -1;
      if (ML_Get_PrintLevel()  >= 10) i = mypid;
      status = ML_CommInfoOP_Deficient_GhostBlk_Check(
                Amatrix->getrow->pre_comm, Amatrix->num_PDEs,i);
      if (status != -1) status = 0;
      else status = 1;
      ML_gsum_scalar_int(&status, &i, comm);
      if ( (mypid == 0) && (status != 0)) {
         printf("**********************************************************\n");
         printf("Switch to the uncoupled aggregation scheme!!!\n\n");
         printf("A deficient ghost block was discovered.  This means that a\n");
         printf("PDE system (numPDEs > 1) is constructed but the ghost part\n");
         printf("of the matrix doesn't conform to the blocking. An example\n");
         printf("with 2x2 blocks follows:\n");
         printf("                   (  x  x       )\n");
         printf("  global view      (  x  x     x )\n");
         printf("  of matrix        (        x  x )\n");
         printf("                   (     x  x  x )\n\n");
         printf("We would have a problem if on two processors this matrix\n");
         printf("is essentially stored as\n");
         printf("  proc 0: (  x  x    )         proc 1: (     x  x )\n");
         printf("          (  x  x  x )                 (  x  x  x )\n\n");
         printf("The problem is that the empty columns have been squeezed\n");
         printf("out and this breaks the 2x2 block structure. Sometimes\n");
         printf("it is possible to trick the system into doing the right\n");
         printf("thing by adding small elements at the right locations.\n");
         printf("For example,\n");
         printf("                   (  x  x  e    )\n");
         printf("  global view      (  x  x     x )\n");
         printf("  of matrix        (  e     x  x )\n");
         printf("                   (     x  x  x )\n\n");
         printf("where e is a small element. This might then force\n");
         printf("something like\n");
         printf("  proc 0: (  x  x  e    )      proc 1: (  e     x  x )\n");
         printf("          (  x  x     x )              (     x  x  x ) .\n\n");
         printf("IT IS IMPORTANT TO RECOGNIZE that this trick might not\n");
         printf("be enough if ghost columns get permuted. In fact, this\n");
         printf("problem can even occur for dense blocks if the order of\n");
         printf("the local ghost columns is different from the order of the\n");
         printf("global columns (ie., consecutive columns within a block\n");
         printf("are no longer consecutive).\n");
         printf("The specific deficient ghost block can be printed by\n");
         printf("setting ML's print level to 10 (or greater).\n");
         printf("**********************************************************\n");
         fflush(stdout);  exit(1);
      }
   }

   /* If the number of points per processor is small, we want to repartition */
   /* repartition the matrix so that only a few processors are used but they */
   /* have more points per processor. This way the aggregation routines will */
   /* work better. The basic idea is to repartition the matrix, update the   */
   /* nullspace to correspond to this new partitioning, compute a tentative  */
   /* prolongator, and permute everything back so that it works with the     */
   /* original A.                                                            */

   switch ( coarsen_scheme )
   {
      case ML_AGGR_UNCOUPLED :
           Ncoarse = ML_Aggregate_CoarsenUncoupled(ag,Amatrix,Pmatrix,comm);
           break;

      case ML_AGGR_COUPLED :
           Ncoarse = ML_Aggregate_CoarsenCoupled(ag,Amatrix,Pmatrix,comm);
           break;

      case ML_AGGR_MIS :
           Ncoarse = ML_Aggregate_CoarsenMIS(ag,Amatrix,Pmatrix,comm);
           break;

      case ML_AGGR_DD :
           Ncoarse = ML_Aggregate_CoarsenDomainDecomp(ag,Amatrix,Pmatrix,comm);
           break;
/*MS*/
      case ML_AGGR_METIS :
           Ncoarse = ML_Aggregate_CoarsenMETIS(ag,Amatrix,Pmatrix,comm);
           break;

      case ML_AGGR_PARMETIS :
           Ncoarse = ML_Aggregate_CoarsenParMETIS(ag,Amatrix,Pmatrix,comm);
           break;

      case ML_AGGR_ZOLTAN :
           Ncoarse = ML_Aggregate_CoarsenZoltan(ag,Amatrix,Pmatrix,comm);
           break;

      case ML_AGGR_USER :
           Ncoarse = ML_Aggregate_CoarsenUser(ag,Amatrix,Pmatrix,comm);
           break;
/*ms*/
/*mgee*/
      case ML_AGGR_VBMETIS :
           Ncoarse = ML_Aggregate_CoarsenVBMETIS(ag,Amatrix,Pmatrix,comm);
           break;

   default :
           if (mypid == 0) printf("ML_Aggregate_Coarsen : invalid scheme.\n");
           exit(1);
           break;
   }


#ifdef ML_DEBUG
   i = 0;
   i = ML_gmax_int(i, comm);
   if ( mypid == 0 && ag->print_flag  < ML_Get_PrintLevel())
      printf("ML_Aggregate_Coarsen ends.\n");
#endif
#ifdef ML_TIMING
   t0 = GetClock() - t0;
   t0 = ML_gsum_double(t0, comm);
   t0 = t0/((double) comm->ML_nprocs);
   if (comm->ML_mypid == 0 && ML_Get_PrintLevel() > 9)
      printf("Aggregation time \t= %e\n",t0);
#endif

   i = -1;
   label = ML_memory_check(NULL);
   if (label != NULL)
     if ( label[0] == 'L')
       if ( ( label[2] == ':') || ( label[3] == ':') )
	 sscanf(&(label[1]),"%d",&i);

   if (ag->keep_agg_information != 1) {
      if ((ag)->aggr_info != NULL) {
	  if ((ag)->aggr_info[ag->cur_level] != NULL)
	    ML_memory_free((void **)&((ag)->aggr_info[ag->cur_level]));
      }
   }
   if (i != 1) ML_memory_check("L%d: agg end",i);
   else ML_memory_check("agg end");

   ag->num_PDE_eqns = Amatrix->num_rigid;

   return Ncoarse;
}

/* ************************************************************************* */
/* Print information about current state of ML_Aggregate                     */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Print( ML_Aggregate *ag )
{
   printf("**************************************************************\n");
   printf("* ML Aggregation information                                 *\n");
   printf("==============================================================\n");
   switch (ag->ordering)
   {
      case 0 : printf("ML_Aggregate : ordering           = natural.\n");
               break;
      case 1 : printf("ML_Aggregate : ordering           = random.\n");
               break;
      case 2 : printf("ML_Aggregate : ordering           = graph.\n");
               break;
   }
   printf("ML_Aggregate : min nodes/aggr     = %d\n",
          ag->min_nodes_per_aggregate);
   printf("ML_Aggregate : max neigh selected = %d\n",
          ag->max_neigh_already_selected);
   switch (ag->attach_scheme)
   {
      case ML_AGGR_MAXLINK :
           printf("ML_Aggregate : attach scheme      = MAXLINK\n");
           break;
      case ML_AGGR_MINRANK :
           printf("ML_Aggregate : attach scheme      = MINRANK\n");
           break;
   }
   /* MS comment this out because it doesn't always work
   switch (ag->coarsen_scheme)
   {
      case ML_AGGR_UNCOUPLED :
           printf("ML_Aggregate : coarsen scheme     = UNCOUPLED\n");
           break;
      case ML_AGGR_COUPLED :
           printf("ML_Aggregate : coarsen scheme     = COUPLED\n");
           break;
   }
   */
   printf("ML_Aggregate : strong threshold   = %e\n", ag->threshold);
   printf("ML_Aggregate : P damping factor   = %e\n",
                                ag->smoothP_damping_factor);
   printf("ML_Aggregate : number of PDEs     = %d\n", ag->num_PDE_eqns);
   printf("ML_Aggregate : number of null vec = %d\n", ag->nullspace_dim);
   printf("ML_Aggregate : smoother drop tol  = %e\n",
          ag->drop_tol_for_smoothing);
   printf("ML_Aggregate : max coarse size    = %d\n", ag->max_coarse_size);
   printf("ML_Aggregate : max no. of levels  = %d\n", ag->max_levels);
   printf("**************************************************************\n");
   return 1;
}

/* ************************************************************************* */
/* Print information about operator complexity                               */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Print_Complexity( ML_Aggregate *ag )
{
   if ( ag->fine_complexity != 0.0 )
      printf("Smoothed Aggregation : operator complexity = %e.\n",
              ag->operator_complexity / ag->fine_complexity);
   else
      printf("Smoothed Aggregation error :  fine complexity = 0.0.\n");
   return 0;
}

/* ************************************************************************* */
/* Exchange boundary function                                                */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_ExchangeBdry(double *vec_data, void *in_comm)
{
   int     N_send_neighbors, msgtype, offset, length;
   int     N_recv_neighbors, fromproc, nbytes;
   int     i, total_send_leng, toproc;
   double  *send_buf = NULL;
   USR_REQ *request;
   ML_Comm *comm;
   ML_Aggregate_Comm *aggr_comm;

   aggr_comm = (ML_Aggregate_Comm *) in_comm;
   N_send_neighbors = aggr_comm->N_send_neighbors;
   N_recv_neighbors = aggr_comm->N_recv_neighbors;
   if (N_send_neighbors == 0 && N_recv_neighbors == 0) return 0;
   comm = aggr_comm->comm;

   nbytes = N_recv_neighbors * sizeof(USR_REQ);
   if (nbytes > 0) ML_memory_alloc( (void **) &request, (unsigned int) nbytes, "AE1" );
   else            request = NULL;
   total_send_leng = 0;
   for ( i = 0; i < N_send_neighbors; i++ )
      total_send_leng += aggr_comm->send_leng[i];
   nbytes = total_send_leng * sizeof(double);
   if (nbytes > 0) ML_memory_alloc( (void **) &send_buf, (unsigned int) nbytes, "AE2" );
   else            send_buf = NULL;
   for ( i = 0; i < total_send_leng; i++ )
   {
      send_buf[i] = vec_data[aggr_comm->send_list[i]];
   }

   /* post receives for all messages */

   offset = aggr_comm->local_nrows;
   for (i = 0; i < N_recv_neighbors; i++)
   {
      msgtype = 1999;
      length  = aggr_comm->recv_leng[i] * sizeof(double);
      fromproc = aggr_comm->recv_neighbors[i];
      comm->USR_irecvbytes((void *) &(vec_data[offset]), (unsigned int) length, &fromproc,
                           &msgtype, comm->USR_comm, request+i);
      offset += aggr_comm->recv_leng[i];
   }

   /* write out all messages */

   offset = 0;
   msgtype = 1999;
   for (i = 0; i < N_send_neighbors; i++)
   {
      length = aggr_comm->send_leng[i] * sizeof(double);
      toproc = aggr_comm->send_neighbors[i];
      comm->USR_sendbytes((void *) &(send_buf[offset]), (unsigned int) length, toproc,
                           msgtype, comm->USR_comm);
      offset += aggr_comm->send_leng[i];
   }

   /* wait for all messages */

   offset = aggr_comm->local_nrows;
   for (i = 0; i < N_recv_neighbors; i++)
   {
      msgtype = 1999;
      length  = aggr_comm->recv_leng[i] * sizeof(double);
      fromproc = aggr_comm->recv_neighbors[i];
      comm->USR_cheapwaitbytes((void *) &(vec_data[offset]), (unsigned int) length, &fromproc,
                           &msgtype, comm->USR_comm, request+i);
      offset += aggr_comm->recv_leng[i];
   }

   ML_memory_free((void**) &request);
   ML_memory_free((void**) &send_buf);
   return 0;
}

/* ************************************************************************* */
/* Exchange data between processors given communication information          */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_ExchangeData(char *recvbuf, char *sendbuf, int N_neighbors,
              int *neighbors, int *recv_leng, int *send_leng, int msgid,
              int datatype, ML_Comm *comm)
{
   int     i, nbytes, fromproc, length, typeleng, msgtype, offset;
   USR_REQ *Request;

   switch ( datatype )
   {
      case ML_CHAR    : typeleng = sizeof(char);   break;
      case ML_INT     : typeleng = sizeof(int);    break;
      case ML_DOUBLE  : typeleng = sizeof(double); break;
      default :         typeleng = datatype;       break;
   }

   nbytes = N_neighbors * sizeof(USR_REQ);
   if ( nbytes > 0 ) ML_memory_alloc( (void **) &Request, (unsigned int) nbytes, "AGZ" );
   else              Request = NULL;
   offset = 0;
   msgtype = msgid;
   for ( i = 0; i < N_neighbors; i++ )
   {
      fromproc = neighbors[i];
      length = recv_leng[i] * typeleng;
      comm->USR_irecvbytes(&recvbuf[offset*typeleng], (unsigned int) length,&fromproc,
#ifdef ML_CPP
                        &msgtype, comm->USR_comm, &Request[i] );
#else
                        &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      offset += recv_leng[i];
   }
   offset = 0;
   msgtype = msgid;
   for ( i = 0; i < N_neighbors; i++ )
   {
      length = send_leng[i] * typeleng;
      comm->USR_sendbytes((void*) &sendbuf[offset*typeleng], (unsigned int) length,
                          neighbors[i], msgtype, comm->USR_comm );
      offset += send_leng[i];
   }
   offset = 0;
   for ( i = 0; i < N_neighbors; i++ )
   {
      fromproc = neighbors[i];
      length = recv_leng[i] * typeleng;
      msgtype = msgid;
      comm->USR_cheapwaitbytes(&recvbuf[offset*typeleng], (unsigned int) length, &fromproc,
#ifdef ML_CPP
                          &msgtype, comm->USR_comm, &Request[i] );
#else
                          &msgtype, comm->USR_comm, (void *) &Request[i] );
#endif
      offset += recv_leng[i];
   }
   if ( Request != NULL ) ML_memory_free((void**) &Request);
   return 0;
}

/* ************************************************************************* */
/* destructor of CSR data structure                                          */
/* ------------------------------------------------------------------------- */

void ML_CSR_MSR_ML_memorydata_Destroy(void *data)
{
   struct ML_CSR_MSRdata *temp;

   temp = (struct ML_CSR_MSRdata *) data;
   if (temp != NULL)
   {
      if (temp->columns != NULL) ML_memory_free( (void**) &(temp->columns));
      if (temp->values  != NULL) ML_memory_free( (void**) &(temp->values));
      if (temp->rowptr  != NULL) ML_memory_free( (void**) &(temp->rowptr));
      ML_memory_free( (void**) &(temp));
   }
}

/* ************************************************************************* */
/* get the aggregation information                                           */
/* ------------------------------------------------------------------------- */

int ML_Gen_Blocks_Aggregates(ML_Aggregate *ag, int level, int *nblocks,
                             int **block_list)
{
   *nblocks = ML_Aggregate_Get_AggrCount( ag, level );
   ML_Aggregate_Get_AggrMap( ag, level, block_list);
   return 0;
}

/* ************************************************************************* */
/* Ray, what does this function do and when is it used.                      */
/*    This looks like something Dawn added. Could we leave it for now        */
/*    until I really figure out what it does?                                */
/* ------------------------------------------------------------------------- */

#ifdef ML_AGGR_PARTEST

int local_to_global_row(int row)
{
   /* global_mapping is an Nrows x 2 array, where the first column
      is the global number of a row and the second column is the local
      number of that row */

   int i, glob_row;

   glob_row=-1;
   i=0;
   while ( i < global_nrows )
   {
      if (global_mapping[i][1] == row)
      {
         glob_row = global_mapping[i][0];
         break;
      }
      i++;
   }
   return(glob_row);
}

#endif

int ML_modified_matvec(void *Amat_in, int ilen, double p[], int olen, double ap[], int num_PDE_eqns)
{
  int i, j;
  int jj, kk;
  double            *p2;
   ML_CommInfoOP     *getrow_comm;
   ML_Operator       *Amat;
   ML_Comm           *comm;
   int allocated_space = 0;
   int *cols = NULL;
   double *vals = NULL;
   int length;
   double *dtemp, diag, best;


   ML_avoid_unused_param((void *) &ilen);

   Amat  = (ML_Operator *) Amat_in;
   comm  = Amat->comm;


   getrow_comm= Amat->getrow->pre_comm;
   if (getrow_comm != NULL) {
      p2 = (double *) ML_allocate((olen+getrow_comm->minimum_vec_size+1)*
                                  sizeof(double));
      dtemp = (double *) ML_allocate((olen+getrow_comm->minimum_vec_size+1)*
                                  sizeof(double));
      for (i = 0; i < olen; i++) p2[i] = p[i];
      ML_exchange_bdry(p2,getrow_comm, olen, comm, ML_OVERWRITE,NULL);
   }
   else {
     p2 = p;
     dtemp = (double *) ML_allocate((olen+1)*sizeof(double));
   }


   for (i = 0; i < olen; i++) {
     ap[i] = -1.0e-20;
     ap[i] = 0.;
     ML_get_matrix_row(Amat, 1, &i , &allocated_space , &cols, &vals,
		       &length, 0);
     ML_random_vec(dtemp,length,Amat->comm);
     diag = -100.;
     best = -100.;
     for (j = 0; j < length; j++) {
       if ((cols[j] == i) && (p2[cols[j]] > 0.)) diag = p2[cols[j]];
       else {
	 if ((dtemp[j] > best) && (p2[cols[j]] > 0.0)) { ap[i] = p2[cols[j]]; best = dtemp[j]; }
       }
       if (diag != -100.) ap[i] = diag;
     }
     /* make sure that all within num_PDE_block are assigned to the same agg */
     if (ap[i] > 0) {
       kk = i/num_PDE_eqns;
       kk = kk*num_PDE_eqns;
       for (jj = kk; jj < kk+num_PDE_eqns; jj++) ap[jj] = ap[i];
     }
   }


  if (getrow_comm != NULL) {
     for (i = 0; i < olen; i++) p[i] = p2[i];
     ML_free(p2);
  }
  ML_free(dtemp);
  if (cols != NULL) ML_free(cols);
  if (vals != NULL) ML_free(vals);
  return(1);
}

int ML_random_global_subset(ML_Operator *Amat, double reduction,
			    int **list, int *length, int num_PDE_eqns)
{
  int Nglobal, itemp;
  int iNtarget, i, *rand_list, *temp_list;
  double Ntarget, dtemp;

  Nglobal = Amat->outvec_leng;
  ML_gsum_scalar_int(&Nglobal,&i,Amat->comm);
  Nglobal = Nglobal/num_PDE_eqns;

  Ntarget = ((double) Nglobal)/reduction;

  Ntarget += .5*( (double) (Nglobal))/(reduction*reduction);
  iNtarget = (int) Ntarget;
      /* add some additional values because we will have */
      /* some duplicate random numbers */

  rand_list = (int *) ML_allocate(iNtarget*(sizeof(int)));
  temp_list = (int *) ML_allocate(iNtarget*(sizeof(int)));

  if (Amat->comm->ML_mypid == 0) {
    for (i = 0; i < iNtarget; i++) {
      ML_random_vec(&dtemp,1,Amat->comm);
      ML_random_vec(&dtemp,1,Amat->comm);
      dtemp = (1 + dtemp);
      rand_list[i] = (int ) (((double) Nglobal)  * dtemp);
      rand_list[i] = (rand_list[i] % Nglobal);
      rand_list[i] = rand_list[i];
      rand_list[i] = rand_list[i];
    }
    ML_az_sort(rand_list, iNtarget, NULL, NULL);
    ML_rm_duplicates(rand_list, &iNtarget);
  }
  else {
    for (i = 0; i < iNtarget; i++) rand_list[i] = 0;
    iNtarget = 0;
  }
  ML_gsum_scalar_int(&iNtarget, &itemp, Amat->comm);
  ML_gsum_vec_int(&rand_list,&temp_list,iNtarget,Amat->comm);
  *length = iNtarget;
  ML_free(temp_list);
  *list = rand_list;

  return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* ML_repartition_matrix() is intended to repartitions a matrix onto a subset*/
/* of processors when the number of rows <= 10 * Nprocs. The basic algorithm */
/* is fairly crude:                                                          */
/*      1) A subset of root unknowns is picked randomly. Each root node      */
/*         will correspond to a processor subdomain.                         */
/*      2) The root nodes are grown into subdomains using a greedy nearest-  */
/*         neighbor algorithm. Repeated multiplies with a modified matrix-   */
/*         vector product routine are used for this. Ties (unknowns adjacent */
/*         to more than one subdomain) are decided randomly.                 */
/*                                                                           */
/* This function is not intended to produce super great partitions. Instead  */
/* its role is to relieve stress from aggregation routines whose quality     */
/* suffers when each processor has roughly one point. The idea is to generate*/
/* P_tent via aggregation on the repartitioned matrix and then permute P_tent*/
/* for use on the original matrix. Specifically, if Q is a permutation       */
/* defining the repartitioning, then                                         */
/*                                                                           */
/*        A_repartition = Q A Q^T                                            */
/* with coarse grid matrix                                                   */
/*        Acoarse_repartition = (P_tent)^T A_repartition P_tent              */
/*                            = (P_tent)^T Q A Q^T P_tent                    */
/* This implies that Q^T P_tent used with the original matrix generates the  */
/* same coarse grid operator.                                                */
/*                                                                           */
/* Note: the time to repartition the matrix is included in its build_time.   */
/*                                                                           */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ML_repartition_matrix(ML_Operator *mat, ML_Operator **new_mat,
              ML_Operator **permutation, ML_Operator **permt,
              int num_PDE_eqns, int Nprocs_ToUse,
              double *xcoord, double *ycoord, double *zcoord,
              int UseImplicitTranspose, ML_Partitioner which_partitioner)
{
 int mypid, nprocs, Nglobal, i, j, the_length;
#if !defined(HAVE_ML_PARMETIS) && !defined(HAVE_ML_ZOLTAN) && !defined(HAVE_ML_JOSTLE)
 int oldj, *the_list = NULL, offset, Nnonzero, oldNnonzero, *itemp = NULL;
  double * d2vec;
#endif
#if defined(ML_TIMING)
  double t0;
#endif

 int *remote_offsets = NULL, *iwork = NULL;
 int *ttemp = NULL, *Procs_WhoGive_TheirUnPermuted_Rows = NULL;
 int Nprocs_WhoGive_TheirUnPermuted_Rows, msgtype, prev;
 int *NRows_IGet_ForPermuted = NULL, *N_UnPermutedRows_ToSend = NULL;
 int Nrows_In_Both_Permuted_And_UnPermuted, fromproc, toproc;
 int Nprocs_WhoGet_MyUnpermutedRows, *Procs_WhoGet_MyUnpermutedRows = NULL;
 int Nrows_Permuted, max_per_processor, MyGlobal_Id, *GlobalId_Rows_ISend=NULL;
 int GlobalId_RowsStaying_WithMe, *columns = NULL, *rowptr = NULL;
 ML_Operator *permt_mat = NULL, *perm_mat = NULL, *permuted_Amat = NULL;
 struct ML_CSR_MSRdata *temp = NULL;
 double *dvec, *values = NULL;
 USR_REQ *request = NULL;
 ML_Comm *comm;
#if defined(HAVE_ML_PARMETIS) || defined(HAVE_ML_ZOLTAN) || defined(HAVE_ML_JOSTLE)
 int *block_list;
#endif

 comm = mat->comm;
 mypid = comm->ML_mypid;
 nprocs = comm->ML_nprocs;
  *new_mat = NULL;
  *permutation = NULL;
  Nglobal = mat->invec_leng;
  ML_gsum_scalar_int(&Nglobal, &i, comm);

#if defined(HAVE_ML_PARMETIS) || defined(HAVE_ML_ZOLTAN) || defined(HAVE_ML_JOSTLE)
   block_list = (int *) ML_allocate(1 + mat->outvec_leng*sizeof(int));
   if (block_list == NULL)
      pr_error("ML_repartition_matrix: out of space\n");

   the_length = Nprocs_ToUse; /* reduces the number of processors with */
                              /* points */

   if (num_PDE_eqns != 1)
     ML_Operator_AmalgamateAndDropWeak(mat, num_PDE_eqns, 0.);
   ML_Operator_BlockPartition(mat, mat->outvec_leng, &the_length,
			      block_list,which_partitioner,xcoord,ycoord,zcoord,num_PDE_eqns);
   if (num_PDE_eqns != 1)
     ML_Operator_UnAmalgamateAndDropWeak(mat, num_PDE_eqns, 0.);

   dvec = (double *) ML_allocate(sizeof(double)*mat->outvec_leng);
   for (i = 0; i < mat->outvec_leng/num_PDE_eqns; i++)
     for (j = 0; j < num_PDE_eqns; j++)
       dvec[j+i*num_PDE_eqns] = (double) ( block_list[i] + 1);

   ML_free(block_list);

#else
   if (mypid == 0 && ML_Get_PrintLevel() > 0)
     printf("ML*WRN* No 3rd party load balancing tool is available.  Continuing with\nML*WRN* repartitioning that uses a simple internal ML method.\n");
  if (Nglobal/num_PDE_eqns >= 12*nprocs) return 1;

  /* Choose a random subset of global unknowns. These will become root nodes */
  /* for a very simple aggregation scheme. These aggregates/subdomains are   */
  /* later assigned to processors defining the new partitioning.             */

  ML_random_global_subset(mat, 5.,&the_list, &the_length, num_PDE_eqns);
  for (i = 0; i < the_length; i++) the_list[i] *= num_PDE_eqns;
#ifdef DEBUG
  if (mypid == 0) {
    for (i = 0; i < the_length; i++)
      printf("%d: ML_reparition_matrix: root_node(%d) = %d\n",mypid,i,the_list[i]);
  }
  fflush(stdout);
#endif

  if (the_length >= nprocs) {
    ML_free(the_list);
    return 1;
  }

  /* Make up a global numbering for the matrix unknowns. The jth  */
  /* unknown on processor k has global id = remote_offsets[k] + j */

  remote_offsets = (int *) ML_allocate(sizeof(int)*(nprocs + 1));
  itemp         = (int *) ML_allocate(sizeof(int)*(Nglobal +
						   nprocs + 1));
  for (i = 0; i < nprocs; i++) remote_offsets[i] = 0;
  remote_offsets[mypid] = mat->invec_leng;
  ML_gsum_vec_int(&remote_offsets,&itemp,nprocs,comm);
  j = 0; oldj = 0;
  for (i = 1; i < nprocs; i++) {
    oldj = j;
    j += remote_offsets[i-1];
    remote_offsets[i-1] = oldj;
  }
  remote_offsets[ nprocs-1] = j;
  remote_offsets[nprocs] = Nglobal;
  offset = remote_offsets[mypid];
#ifdef DEBUG
  if (mypid == 0)
    for (i = 0; i <= nprocs; i++)
      printf("ML_repartition_matrix: REMOTE OFFSET(%d) = %d\n",i,remote_offsets[i]);
  fflush(stdout);
#endif

  /* Create a vector which initially has all elements set */
  /* to zero except those corresponding to the root nodes */
  /* defined in 'the_list'. These are set to the position */
  /* in 'the_list' where they reside. A modified matrix-  */
  /* vector product routine follows using this vector.    */
  /* This routine will effectively grow the root nodes    */
  /* into domains.                                        */

  /* Make sure there is enough space for ghost nodes.     */
  i = mat->invec_leng + 1;
  if (mat->getrow->pre_comm != NULL)
    i += mat->getrow->pre_comm->minimum_vec_size;

  dvec = (double *) ML_allocate(sizeof(double)*i);
  d2vec = (double *) ML_allocate(sizeof(double)*i);

  for (i = 0; i < mat->invec_leng; i++) dvec[i] = 0.;
  Nnonzero = the_length;
  for (i = 0; i < the_length; i++) {
    if ( (the_list[i] >= offset) &&
	 (the_list[i] < offset + mat->invec_leng)){
      for (j = 0 ; j < num_PDE_eqns; j++)
	dvec[the_list[i]-offset+j] = (double ) (i + 1);
    }
  }

  /* Grow the subdomains by using a modified matvec routine. This */
  /* routine takes dvec and produces d2vec. This routine works as */
  /* follows:                                                     */
  /*    if dvec[k] != 0              ========> d2vec[k] = dvec[k] */
  /*    if dvec[k] == 0 and                                       */
  /*	  A(k,j) != 0 for some                                    */
  /*      dvec[j]!= 0                ========> d2vec[k] = dvec[j] */
  /*    else                         ========> d2vec[k] = 0       */
  /* NOTE: if several A(k,j)'s are nonzero with dvec[j] nonzero,  */
  /* then one j is chosen at random.                              */

  oldNnonzero = 0;
  while ( (Nglobal - Nnonzero > 0) && (Nnonzero > oldNnonzero)) {
    ML_modified_matvec( mat, mat->invec_leng, dvec,
			mat->outvec_leng, d2vec, num_PDE_eqns);
    oldNnonzero = Nnonzero;
    Nnonzero = 0;
    for (i = 0; i < mat->outvec_leng; i++) {
      dvec[i] = d2vec[i];
      if (dvec[i] > 0.0) Nnonzero++;
    }
    ML_gsum_scalar_int(&Nnonzero,&i,comm);
  }
  if (d2vec != NULL) ML_free(d2vec);

  /* Assign any singletons to a new domain */

  for (i = 0; i < mat->outvec_leng; i++) {
    if ( dvec[i] == 0.) dvec[i] = (double) (the_length + 1);
  }

#ifdef DEBUG
  for (i = 0; i < mat->outvec_leng; i++)
    printf("%d: ML_repartition_matrix: domain_assignment(%d) = %e\n",mypid,i,dvec[i]);
  fflush(stdout);
#endif

#endif
  /* Change dvec so that it now contains the global id (instead */
  /* of the subdomain number) corresponding to the row in the   */
  /* permuted system.                                           */

  /* Compute the number of processors (not including myself) */
  /* owning rows in the unpermuted system that correspond    */
  /* to rows that I will now own in the permuted system.     */

  iwork  = (int *) ML_allocate(sizeof(int)*nprocs);
  ttemp  = (int *) ML_allocate(sizeof(int)*nprocs);
  for (i = 0; i < nprocs; i++) iwork[i] = 0;
  for (i = 0; i < mat->invec_leng; i++) {
    j = (int) dvec[i];
    j--;
    if (j < 0) j = 0;
    if ( j > nprocs)
      pr_error("%d: dvec > Nprocs?  %d\n",mypid,(int) dvec[i]);
    iwork[j] = 1;
  }
  iwork[mypid] = 0;
  ML_gsum_vec_int(&iwork,&ttemp,nprocs,comm);
  Nprocs_WhoGive_TheirUnPermuted_Rows = iwork[mypid];
  if (ttemp != NULL) ML_free(ttemp);
  if (iwork != NULL) ML_free(iwork);


  /* Compute the # of my unpermuted rows that will be   */
  /* assigned to each processor in the permuted system. */

  if (Nprocs_WhoGive_TheirUnPermuted_Rows > 0)
    request = (USR_REQ *) ML_allocate(Nprocs_WhoGive_TheirUnPermuted_Rows*
				      sizeof(USR_REQ));
  else            request = NULL;
  Procs_WhoGive_TheirUnPermuted_Rows=(int *) ML_allocate(sizeof(int)*
							 (Nprocs_WhoGive_TheirUnPermuted_Rows+1));
  NRows_IGet_ForPermuted = (int *) ML_allocate(sizeof(int)*
					       (Nprocs_WhoGive_TheirUnPermuted_Rows+1));
  N_UnPermutedRows_ToSend = (int *) ML_allocate(sizeof(int)*(nprocs));

  for (i = 0; i < nprocs; i++) N_UnPermutedRows_ToSend[i] = 0;
  for (i = 0; i < mat->invec_leng; i++) {
    j = (int) dvec[i];
    j--;
    if (j < 0) j = 0;
    (N_UnPermutedRows_ToSend[j])++;
  }
  Nrows_In_Both_Permuted_And_UnPermuted = N_UnPermutedRows_ToSend[mypid];

  /* Now send the # of my unpermuted rows that will be    */
  /* assigned to other processors in the permuted system. */

  for (i = 0; i < Nprocs_WhoGive_TheirUnPermuted_Rows; i++)    {
    msgtype = 1794;
    fromproc = -1;
    comm->USR_irecvbytes((void *)&(NRows_IGet_ForPermuted[i]),
			 sizeof(int),&fromproc,&msgtype,comm->USR_comm,request+i);
    if (fromproc != -1) Procs_WhoGive_TheirUnPermuted_Rows[i] = fromproc;
  }
  Nprocs_WhoGet_MyUnpermutedRows = 0;
  for (i = 0; i < nprocs; i++) {
    if ( (N_UnPermutedRows_ToSend[i] != 0) && (i != mypid)) {
      msgtype = 1794;
      toproc = i;
      comm->USR_sendbytes((void *) &(N_UnPermutedRows_ToSend[i]),sizeof(int),
			  toproc,msgtype, comm->USR_comm);
      Nprocs_WhoGet_MyUnpermutedRows++;
    }
  }
  Procs_WhoGet_MyUnpermutedRows = (int *) ML_allocate(sizeof(int)*
						      (Nprocs_WhoGet_MyUnpermutedRows+1));
  Nprocs_WhoGet_MyUnpermutedRows = 0;
  for (i = 0; i < nprocs; i++) {
    if ( (N_UnPermutedRows_ToSend[i] != 0) && (i != mypid))
      Procs_WhoGet_MyUnpermutedRows[Nprocs_WhoGet_MyUnpermutedRows++] = i;
  }
  if (N_UnPermutedRows_ToSend != NULL) ML_free(N_UnPermutedRows_ToSend);

  for (i = 0; i < Nprocs_WhoGive_TheirUnPermuted_Rows; i++)    {
    msgtype = 1794;
    fromproc = -1;
    comm->USR_cheapwaitbytes((void *) &(NRows_IGet_ForPermuted[i]),sizeof(int),
			     &fromproc,&msgtype, comm->USR_comm, request+i);
    if (fromproc != -1) Procs_WhoGive_TheirUnPermuted_Rows[i] = fromproc;
  }
  ML_az_sort(Procs_WhoGive_TheirUnPermuted_Rows,
	     Nprocs_WhoGive_TheirUnPermuted_Rows,NRows_IGet_ForPermuted,NULL);

  Nrows_Permuted = 0;
  for (i = 0; i < Nprocs_WhoGive_TheirUnPermuted_Rows; i++)
    Nrows_Permuted += NRows_IGet_ForPermuted[i];
  Nrows_Permuted += Nrows_In_Both_Permuted_And_UnPermuted;
  max_per_processor = ML_gmax_int(Nrows_Permuted, comm);

  MyGlobal_Id = max_per_processor*mypid;

  if (Nprocs_WhoGet_MyUnpermutedRows > Nprocs_WhoGive_TheirUnPermuted_Rows ) {
    if (request != NULL) ML_free(request);
    request = (USR_REQ *) ML_allocate(Nprocs_WhoGet_MyUnpermutedRows*
				      sizeof(USR_REQ));
  }

  /* send back the start index for each of rows represented in dvec */
  /* so that we can convert dvec to a global id                     */
  GlobalId_Rows_ISend = (int *) ML_allocate(sizeof(int)*(nprocs));
  for (i = 0; i < nprocs; i++) GlobalId_Rows_ISend[i] = -1;


  for (i = 0; i < Nprocs_WhoGet_MyUnpermutedRows; i++)    {
    msgtype = 4971;
    fromproc = Procs_WhoGet_MyUnpermutedRows[i];
    comm->USR_irecvbytes((void *) &(GlobalId_Rows_ISend[fromproc]),sizeof(int),
			 &fromproc,&msgtype, comm->USR_comm, request+i);
  }
  prev = -1;
  GlobalId_RowsStaying_WithMe = MyGlobal_Id;

  for (i = 0; i < Nprocs_WhoGive_TheirUnPermuted_Rows; i++)    {
    msgtype = 4971;
    toproc = Procs_WhoGive_TheirUnPermuted_Rows[i];

    if ( (prev < mypid) && (toproc > mypid) ) {
      MyGlobal_Id += Nrows_In_Both_Permuted_And_UnPermuted;
    }
    prev = toproc;
    comm->USR_sendbytes((void *) &MyGlobal_Id, sizeof(int), toproc,
			msgtype, comm->USR_comm);
    MyGlobal_Id += NRows_IGet_ForPermuted[i];
    if (toproc < mypid) {
      GlobalId_RowsStaying_WithMe = MyGlobal_Id;
    }
  }
  if (NRows_IGet_ForPermuted != NULL) ML_free(NRows_IGet_ForPermuted);
  if (Procs_WhoGive_TheirUnPermuted_Rows != NULL)
    ML_free(Procs_WhoGive_TheirUnPermuted_Rows);

  for (i = 0; i < Nprocs_WhoGet_MyUnpermutedRows; i++)    {
    msgtype = 4971;
    fromproc = Procs_WhoGet_MyUnpermutedRows[i];
    comm->USR_cheapwaitbytes((void *) &(GlobalId_Rows_ISend[fromproc]),
			     sizeof(int),&fromproc,&msgtype,comm->USR_comm,request+i);
  }
  if (Procs_WhoGet_MyUnpermutedRows != NULL)
    ML_free(Procs_WhoGet_MyUnpermutedRows);
  if (request != NULL) ML_free(request);

  for (i = 0; i < mat->invec_leng; i++) {
    j = (int) dvec[i];
    j--;
    if (j < 0) j = 0;
    if (j != mypid) {
      dvec[i] = GlobalId_Rows_ISend[j];
      GlobalId_Rows_ISend[j]++;
    }
    else dvec[i] = GlobalId_RowsStaying_WithMe++;
  }

  if (GlobalId_Rows_ISend != NULL) ML_free(GlobalId_Rows_ISend);
  /* build the new matrix */
  columns = (int *) ML_allocate(sizeof(int)*(mat->invec_leng+1));
  values  = (double *) ML_allocate(sizeof(double)*(mat->invec_leng+1));
  rowptr = (int *) ML_allocate(sizeof(int)*(mat->invec_leng+1));
  rowptr[0] = 0;
  for (i = 0; i < mat->invec_leng; i++) {
    rowptr[i+1] = i+1;
    values[i]   = 1.;
    columns[i]  = (int) dvec[i];
  }
  temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  temp->columns = columns;
  temp->values  = values;
  temp->rowptr  = rowptr;
  permt_mat = ML_Operator_Create(comm);
  ML_Operator_Set_ApplyFuncData(permt_mat,Nrows_Permuted,mat->invec_leng,
				temp,Nrows_Permuted,NULL,0);
  ML_Operator_Set_Getrow(permt_mat, mat->invec_leng, CSR_getrow);
  ML_Operator_Set_ApplyFunc( permt_mat, CSR_matvec);
  permt_mat->data_destroy = ML_CSR_MSRdata_Destroy;
  permt_mat->max_nz_per_row = 1;
  permt_mat->min_nz_per_row = 1;
  permt_mat->N_nonzeros     = mat->invec_leng;
  *permt = ML_Operator_Create(comm);
  ML_back_to_csrlocal(permt_mat, *permt, max_per_processor);
  ML_Operator_Destroy(&permt_mat);
  ML_Operator_ChangeToChar(*permt);

  /* Check if this matrix is low storage character format */
  /* is used. In this case only, it is possible to get    */
  /* rid of rowptr as there is an option to not use it    */

  if ( ((*permt)->getrow->func_ptr == cCSR_getrows) &&
       ((*permt)->matvec->func_ptr == cCSR_matvec)) {
    temp = (struct ML_CSR_MSRdata *) (*permt)->data;
    ML_free(temp->rowptr);
    temp->rowptr = NULL;
  }

  perm_mat = ML_Operator_Create(comm);
  ML_Operator_Set_Label(perm_mat,"Rebalancing matrix");
  ML_Operator_Transpose_byrow(*permt, perm_mat);
  ML_Operator_ChangeToChar(perm_mat);

  /* Check if this matrix is low storage character format */
  /* is used. In this case only, it is possible to get    */
  /* rid of rowptr as there is an option to not use it    */

  if ( ((perm_mat)->getrow->func_ptr == cCSR_getrows) &&
       ((perm_mat)->matvec->func_ptr == cCSR_matvec)) {
    temp = (struct ML_CSR_MSRdata *) (perm_mat)->data;
    ML_free(temp->rowptr);
    temp->rowptr = NULL;
  }

  if (ML_Get_PrintLevel() > 15) {
    printf("(pid %d, level %d): permutation matrix has %10d rows, %10d columns\n",mypid,mat->to->levelnum,perm_mat->outvec_leng,
	       perm_mat->invec_leng);
  }

  permuted_Amat = ML_Operator_Create(comm);

#if defined(ML_TIMING)
  t0 = GetClock();
#endif
  ML_rap(perm_mat,mat,*permt,permuted_Amat, ML_CSR_MATRIX);
#if defined(ML_TIMING)
  t0 = GetClock() - t0;
#endif
  ML_Operator_Copy_Statistics(mat,permuted_Amat);
#if defined(ML_TIMING)
  permuted_Amat->build_time += t0;
#endif
  if (UseImplicitTranspose)
    ML_Operator_ImplicitTranspose(perm_mat,*permt, ML_FALSE);

  *new_mat = permuted_Amat;
  *permutation = perm_mat;
  ML_Operator_Set_Label(*permutation,"perm");
  ML_Operator_Set_Label(*permt,"permt");

  if (dvec != NULL) ML_free(dvec);
  if (remote_offsets != NULL) ML_free(remote_offsets);

  return 0;
} /*ML_repartition_matrix()*/

/* ************************************************************************* */
/* Use A instead of AT in restriction                                        */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_SmoothRestrictionWithA( ML_Aggregate *ag )
{
  ag->Restriction_smoothagg_transpose = ML_TRUE;
  return 0;
}

/* ************************************************************************* */
/* Use A in restriction (default)                                            */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_SmoothRestrictionWithAT( ML_Aggregate *ag )
{
  ag->Restriction_smoothagg_transpose = ML_FALSE;
  return 0;
}

/* ************************************************************************* */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_NodalCoordinates(ML* ml, ML_Aggregate *ag, double *ptr)
{
  int i;
  int MaxLevels = ml->ML_num_levels;
  assert (MaxLevels);
  assert (ptr != 0);

  if (ag->nodal_coord)
    ML_free(ag->nodal_coord);

  ag->nodal_coord = (double**) ML_allocate((sizeof(double*) * (MaxLevels)));
  assert (ag->nodal_coord != NULL);
  for (i = 0 ; i < MaxLevels ; ++i)
    ag->nodal_coord[i] = NULL;

  /* NOTE: fine-grid is ALWAYS at position 0.
   * ML_DECREASING can be handle by passing relative level number */
  ag->nodal_coord[0] = ptr;

  return 0;
}

/* ************************************************************************* */
/* Set the number of dimensions                                              */
/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_Dimensions(ML_Aggregate *ag, int N_dimensions)
{
  assert (N_dimensions > 0);
  assert (N_dimensions < 4);
  ag->N_dimensions = N_dimensions;
  return 0;
}


/*******************************************************************************
 Function ML_repartition_Acoarse

    Summary: Redistributes a matrix across a (sub)set of processors.

    Input:
      ReturnPerm  int     if ML_TRUE (1), destroy permutation matrices
                          if ML_FALSE (0), return permutation matrices
                            (it is caller's responsibility to destroy them)

    Output:
      NULL                , if ReturnPerm == ML_FALSE
      array of *ML_Operator & length 2, if ReturnPerm == ML_TRUE
*******************************************************************************/

ML_Operator** ML_repartition_Acoarse(ML *ml, int fine, int coarse,
               ML_Aggregate *ag, int R_is_Ptranspose, int ReturnPerm)
{
  ML_Operator *Amatrix, *Rmat, *Pmat, *perm, *permt, *newA, *newP, *newR;
  ML_Operator **permvec=NULL;
  int status, offset1, offset2, j, flag = 0;
  double *new_null;
  int ml_gmin, ml_gmax, Nprocs_ToUse;
  double ml_gsum;
  double *xcoord = NULL, *ycoord = NULL, *zcoord = NULL;
  double *new_xcoord, *new_ycoord, *new_zcoord;
  int UseImplicitTranspose;
  ML_Partitioner which_partitioner;
  ML_Aggregate_Viz_Stats *grid_info;
  int N_dimensions=0;
  int haveCoordinates = 0;
  double t0=0.0,delta=0.0;

  StartTimer(&t0);


  if (ML_Repartition_Status(ml) == ML_FALSE)
    return NULL;

  Amatrix = &(ml->Amat[coarse]);
  Rmat = &(ml->Rmat[fine]);
  Pmat = &(ml->Pmat[coarse]);

  if ((ml->MinPerProc_repartition == -1) &&
      (ml->LargestMinMaxRatio_repartition == -1.) &&
      (ml->PutOnSingleProc_repartition == -1.))
    return NULL;

  /*FIXME (JJH) for a rectangular matrix, invec_leng != outvec_leng .... */
  ml_gmax = ML_gmax_int(Amatrix->invec_leng,ml->comm);
  ml_gmin = Amatrix->invec_leng;
  if (ml_gmin == 0) ml_gmin = ml_gmax; /* don't count */
                                      /* empty processors */
  ml_gmin = ML_gmin_int(ml_gmin,ml->comm);
  ml_gsum = ML_gsum_double((double)Amatrix->invec_leng,ml->comm);

  if ( (ml->MinPerProc_repartition != -1) &&
       (ml_gmin < ml->MinPerProc_repartition))
    flag = 1;
  else if ((ml->LargestMinMaxRatio_repartition != -1.) &&
         ((((double) ml_gmax)/((double) ml_gmin)) >
          ml->LargestMinMaxRatio_repartition))
    flag = 1;


  Nprocs_ToUse = ml->comm->ML_nprocs;
  if ( (flag == 1) && (ml->MinPerProc_repartition != -1))
  {
    /* compute how many processors to use in the repartitioning */
    double ttt = ml_gsum/ml->MinPerProc_repartition;
    if (ttt > ml->comm->ML_nprocs)
      ttt = ml->comm->ML_nprocs;
    if (ttt < 1) ttt = 1;
    Nprocs_ToUse = (int) ttt;
  }

  if (ml_gsum < ml->PutOnSingleProc_repartition) {
    flag = 1;
    Nprocs_ToUse = 1;
  }

  if (flag == 0)
    return NULL;

  if (ML_Get_PrintLevel() > 0 && !ml->comm->ML_mypid) {
    printf("Repartitioning (level %d): min rows per proc = %d\n",coarse,ml->MinPerProc_repartition);
    printf("Repartitioning (level %d): largest max/min ratio = %2.3e\n",coarse,ml->LargestMinMaxRatio_repartition);
    printf("Repartitioning (level %d): max #rows (global) that fits on one proc = %d\n",coarse,ml->PutOnSingleProc_repartition);
    printf("Repartitioning (level %d): #proc to use in repartitioning = %d\n",coarse,Nprocs_ToUse);
  }

  grid_info = (ML_Aggregate_Viz_Stats *) ml->Grid[coarse].Grid;
  which_partitioner = ML_Repartition_Get_Partitioner(ml);
  /* no coordinates supplied */
  if (grid_info == NULL) {
    haveCoordinates = 0;
    if ((ml->comm->ML_mypid == 0)
        && (ML_Get_PrintLevel() > 0)
        && (which_partitioner == ML_USEZOLTAN)
        && (Nprocs_ToUse > 1))
      printf("ML*WRN* No grid structure found. This is not necessarily an\nML*WRN* error, but repartitioning with Zoltan is impossible.\n\n");
  }
  else if (grid_info->x == NULL || grid_info->y == NULL) {
    haveCoordinates = 0;
    if ((ml->comm->ML_mypid == 0)
        && (ML_Get_PrintLevel() > 0)
        && (which_partitioner == ML_USEZOLTAN)
        && (Nprocs_ToUse > 1))
      printf("ML*WRN* Either x- or y-coordinates are missing. This is not necessarily an\nML*WRN* error, but repartitioning with Zoltan is impossible.\n\n");
  }
  else {
    haveCoordinates = 1;
    xcoord = grid_info->x;
    ycoord = grid_info->y;
    zcoord = grid_info->z;
    N_dimensions = grid_info->Ndim;
    if (N_dimensions < 1 || N_dimensions > 3) {
      N_dimensions = 0;
      if (xcoord != NULL) N_dimensions++;
      if (ycoord != NULL) N_dimensions++;
      if (zcoord != NULL) N_dimensions++;
      grid_info->Ndim = N_dimensions;
      if ((ml->comm->ML_mypid == 0) && (ML_Get_PrintLevel() > 0))
        printf("ML*WRN* ML_repartition_Acoarse: problem dimension was not previously set.\nML*WRN* Now setting dimension to %d.\n",N_dimensions);
    }
    if(ag != NULL){
      if (ag->N_dimensions != N_dimensions) {
	N_dimensions = ag->N_dimensions;
	if  (N_dimensions < 3)  zcoord = NULL;
	if  (N_dimensions < 2)  ycoord = NULL;
      }
    }
  }

  /* Turn off implicit transpose because the getrow is needed to apply
     the permutation matrices. */
  if (ReturnPerm == ML_TRUE) UseImplicitTranspose = ML_FALSE;
  else UseImplicitTranspose = ML_TRUE;

  status = ML_repartition_matrix(Amatrix, &newA, &perm, &permt,
                 Amatrix->num_PDEs, Nprocs_ToUse,
                 xcoord, ycoord, zcoord, UseImplicitTranspose,
                                 which_partitioner);

  if (status == 0)
  {
    int i;
    double * tmp_coord = 0;
    if (ag !=NULL)
    {
     if (ag->nullspace_vect != NULL) {
       new_null = (double *) ML_allocate(sizeof(double)*ag->nullspace_dim*
                  perm->outvec_leng);
       offset1 = 0;
       offset2 = 0;
       for (j = 0; j < ag->nullspace_dim; j++)
       {
         ML_Operator_Apply(perm, perm->invec_leng,
               &((ag->nullspace_vect)[offset1]),
               perm->outvec_leng, &(new_null[offset2]));

         offset1 += perm->invec_leng;
         offset2 += perm->outvec_leng;
       }
       ML_Aggregate_Set_NullSpace(ag, ag->num_PDE_eqns, ag->nullspace_dim,
                                  new_null,perm->outvec_leng);
       ML_free(new_null);
     } /* if (ag->nullspace_vect != NULL) */
    } /* if (ag !=NULL) */

    tmp_coord = (double *) ML_allocate(sizeof(double)*(perm->invec_leng +1));
    if (xcoord != NULL) {
      new_xcoord = (double *) ML_allocate(sizeof(double)*(N_dimensions)*
                                        (perm->outvec_leng +1));

      /* make sure coordinate is setup as a degree of freedom vector */
      for (i=0; i < perm->invec_leng ; i++) {
         tmp_coord[i] = xcoord[i/Amatrix->num_PDEs];
      }

      ML_Operator_Apply(perm, perm->invec_leng,
            tmp_coord, perm->outvec_leng, new_xcoord);
      ML_free(grid_info->x);
      grid_info->x = new_xcoord;
    }
    if (ycoord != NULL) {
      new_ycoord = (double *) ML_allocate(sizeof(double)*(N_dimensions)*
                                          (perm->outvec_leng +1));

      /* make sure coordinate is setup as a degree of freedom vector */
      for (i=0; i < perm->invec_leng ; i++) {
         tmp_coord[i] = ycoord[i/Amatrix->num_PDEs];
      }

      ML_Operator_Apply(perm,
                        perm->invec_leng, tmp_coord,
                        perm->outvec_leng, new_ycoord);
      ML_free(grid_info->y);
      grid_info->y = new_ycoord;
    }
    if (zcoord != NULL) {
      new_zcoord = (double *) ML_allocate(sizeof(double)*(N_dimensions)*
                                          (perm->outvec_leng +1));

      /* make sure coordinate is setup as a degree of freedom vector */
      for (i=0; i < perm->invec_leng ; i++) {
         tmp_coord[i] = zcoord[i/Amatrix->num_PDEs];
      }

      ML_Operator_Apply(perm,
                        perm->invec_leng, tmp_coord,
                        perm->outvec_leng, new_zcoord);
      ML_free(grid_info->z);
      grid_info->z = new_zcoord;
    }
    ML_free(tmp_coord);

    ML_Operator_Move2HierarchyAndDestroy(&newA, Amatrix);

    /* start of MS modif, 26-Mar-05                         */
    /* This lines are required by ML_Project_Coordinates(), */
    /* in ML_Gen_MultiLevelHierarchy(), since I use Ptent   */
    /* to project the coordinates down to the next level    */
    if (ag != NULL) {
      if (ag->P_tentative != 0) {
        if (ag->P_tentative[coarse] != 0)
        {
          newP = ML_Operator_Create(Pmat->comm);
          ML_Operator_Set_Label(ag->P_tentative[coarse],"rebalance Ptentative");
          ML_2matmult(ag->P_tentative[coarse], permt, newP, ML_CSR_MATRIX);
          ML_Operator_Destroy(&(ag->P_tentative[coarse]));
          ag->P_tentative[coarse] = newP;
          newP = NULL;
        }
      }
    }
    /* end of MS modif */

    /* do a mat-mat mult to get the appropriate P for the */
    /* unpartitioned matrix.                              */

    newP = ML_Operator_Create(Pmat->comm);
    ML_2matmult(Pmat, permt, newP, ML_CSR_MATRIX);
    ML_Operator_Copy_Statistics(Pmat,newP);
    ML_Operator_Move2HierarchyAndDestroy(&newP, Pmat);

    /* used to print error message for those         */
    /* attempting to semicoarsen afer repartiitoning.*/
    if (ml->ML_finest_level == 0) {
      for (i=coarse; i < ml->ML_num_levels; i++) {(ml->Pmat)[i].NumZDir = -7;(ml->Pmat)[i].Zorientation= 1;}
    }
    else {
      for (i=coarse; i >= 0; i--) {(ml->Pmat)[i].NumZDir = -7;(ml->Pmat)[i].Zorientation= 1;}
    }


    if (R_is_Ptranspose == ML_TRUE) {
      newR = ML_Operator_Create(Rmat->comm);
      ML_Operator_Transpose(Pmat, newR);
      ML_Operator_Move2HierarchyAndDestroy(&newR, Rmat);
    }
    else if (Rmat->getrow->post_comm == NULL) {
      newR = ML_Operator_Create(Rmat->comm);
      ML_Operator_Set_Label(perm,"rebalance R");
      ML_2matmult(perm, Rmat, newR, ML_CSR_MATRIX);
      ML_Operator_Copy_Statistics(Rmat,newR);
      ML_Operator_Move2HierarchyAndDestroy(&newR, Rmat);
    }
    else {
      printf("ML_repartition_Acoarse: 2matmult does not work properly if\n");
      printf("   rightmost matrix in multiply is created with an implicit\n");
      printf("   transpose (e.g. ML_Gen_Restrictor_TransP). If R is P^T,\n");
      printf("   then invoke as ML_repartition_Acoarse(..., ML_TRUE). If\n");
      printf("   R is not P^T but an implicit transpose is used, then try\n");
      printf("   to remove implicit transpose with: \n\n");
      printf("   ML_Operator_Transpose_byrow( &(ml->Pmat[next]),&(ml->Rmat[level]));\n");
      printf("   ML_Operator_Set_1Levels(&(ml->Rmat[level]),&(ml->SingleLevel[level]), &(ml->SingleLevel[next]));\n");
      exit(1);
    }
    if (ReturnPerm == ML_FALSE) {
      ML_Operator_Destroy(&perm);
      ML_Operator_Destroy(&permt);
    }
    else {
      permvec = (ML_Operator **) ML_allocate(2 * sizeof(ML_Operator*));
      ML_allocate_check(permvec);
      permvec[0] = perm;
      permvec[1] = permt;
    }
    ML_Operator_ChangeToSinglePrecision(&(ml->Pmat[coarse]));
  } /*if (status == 0)*/

  StopTimer(&t0,&delta);
  if (ML_Get_PrintLevel() > 9)
    ReportTimer(delta,"Time spent in ML_repartition_Acoarse",ml->comm);

  if (ReturnPerm == ML_FALSE)
    return NULL;
  else
    return permvec;
} /*ML_repartition_Acoarse()*/
