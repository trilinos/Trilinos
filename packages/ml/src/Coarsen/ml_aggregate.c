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
#define ML_AGGR_UNCOUPLED    1
#define ML_AGGR_COUPLED      2
#define ML_AGGR_MIS          3
#define ML_AGGR_DD           4
#define ML_AGGR_HYBRIDUC     5
#define ML_AGGR_HYBRIDUM     6
/*MS*/
#define ML_AGGR_METIS        7
#define ML_AGGR_PARMETIS     8
/*ms*/

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
   (*ag)->spectral_radius_scheme     = 1;  /* compute it */
   (*ag)->smoothP_type               = 0;  /* point type */
   (*ag)->num_PDE_eqns               = 1;
   (*ag)->nullspace_dim              = 1;
   (*ag)->nullspace_vect             = NULL;
   (*ag)->nullspace_corrupted        = ML_EMPTY;
   (*ag)->max_levels                 = 0;
   (*ag)->max_coarse_size            = 32;
   (*ag)->begin_level                = 0;
   (*ag)->cur_level                  = 0;
   (*ag)->aggr_info                  = NULL;
   (*ag)->aggr_count                 = NULL;
   (*ag)->drop_tol_for_smoothing     = 0.0;
   (*ag)->fine_complexity            = 0.0;
   (*ag)->nvblocks                   = 0;
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
   (*ag)->field_of_values            = NULL;
/*MS*/
   (*ag)->block_scaled_SA            = 0;

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

      if( (*ag)->field_of_values != NULL ) ML_free( (*ag)->field_of_values );
/*MS*/
      ML_memory_free( (void **) ag );
      (*ag) = NULL;

   }
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

int ML_Aggregate_Set_CoarsenSchemeLevel_METIS( int level, int MaxLevels,
					       ML_Aggregate *ag  )
{
  return( ML_Aggregate_Set_CoarsenSchemeLevel(level, MaxLevels,
					      ag, ML_AGGR_METIS) );
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_CoarsenSchemeLevel_ParMETIS( int level, int MaxLevels,
						  ML_Aggregate *ag  )
{
  return( ML_Aggregate_Set_CoarsenSchemeLevel(level, MaxLevels,
					      ag, ML_AGGR_PARMETIS) );
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
/* set method to estimate spectral radius of A                          */
/* -------------------------------------------------------------------- */

int ML_Aggregate_Set_SpectralNormScheme_Calc( ML_Aggregate *ag )
{
   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      printf("ML_Aggregate_Set_SpectralNormScheme_Calc : wrong object. \n");
      exit(-1);
   }
   ag->spectral_radius_scheme = 1;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_SpectralNormScheme_Anorm( ML_Aggregate *ag )
{
   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      printf("ML_Aggregate_Set_SpectralNormScheme_Anorm : wrong object. \n");
      exit(-1);
   }
   ag->spectral_radius_scheme = 0;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_SpectralNormScheme_Anasazi( ML_Aggregate *ag)
{
   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      printf("ML_Aggregate_Set_SpectralNormScheme_Anorm : wrong object. \n");
      exit(-1);
   }
   ag->spectral_radius_scheme = 2;
   return 0;
}

/* ------------------------------------------------------------------------- */

int ML_Aggregate_Set_SpectralNormScheme_PowerMethod( ML_Aggregate *ag)
{
   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      printf("ML_Aggregate_Set_SpectralNormScheme_PowerMethod : wrong object. \n");
      exit(-1);
   }
   ag->spectral_radius_scheme = 3;
   return 0;
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
	
      ML_memory_alloc((void **)&(ag->nullspace_vect), nbytes, "ns");

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
	
     ML_memory_alloc((void **)&(ag->nullspace_vect), nbytes, "ns");
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
   int i=1, ndofs, Ncoarse, coarsen_scheme, j, offset1, offset2, status;
   int mypid, nprocs;
   double *new_null;
   ML_Operator *newA = NULL, *permt = NULL, *perm = NULL, *oldA = NULL;
   ML_Operator *newP = NULL, *oldP = NULL;

#ifdef ML_TIMING
   double t0;
   t0 = GetClock();
#endif

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
   if ( ndofs < 2 ) ndofs = 0; else ndofs = 1;
   ML_gsum_scalar_int(&ndofs, &i, comm);
   /*
   gmin = ML_gmax_double((double) (-1.0 * Amatrix->outvec_leng) , comm);
   if (comm->ML_mypid == 0)
      printf("Smallest Amatrix->outvec_leng = %d\n", (int) (-1 * gmin));
   */
   nprocs = comm->ML_nprocs;
   /*MS*/
   if( ag->coarsen_scheme_level == NULL ) {
     coarsen_scheme = ag->coarsen_scheme;
   } else {
     coarsen_scheme = ag->coarsen_scheme_level[ag->cur_level];
   }
   /*ms*/
   
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
/*ms*/
      else 
         coarsen_scheme = ML_AGGR_UNCOUPLED;
   }
   else 
   {
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
/*ms*/
      else
      {
         (*Pmatrix) = NULL;
         return 0;
      }
   }

   /* If the number of points per processor is small, we want to repartition */
   /* repartition the matrix so that only a few processors are used but they */
   /* have more points per processor. This way the aggregation routines will */
   /* work better. The basic idea is to repartition the matrix, update the   */
   /* nullspace to correspond to this new partitioning, compute a tentative  */
   /* prolongator, and permute everything back so that it works with the     */
   /* original A.                                                            */

#ifdef ML_repartition
   status = ML_repartition_matrix(Amatrix, &newA, &perm, &permt, ag->num_PDE_eqns);
   
   if (status == 0) {
     if (ag->nullspace_vect != NULL) {
       new_null = ML_allocate(sizeof(double)*ag->nullspace_dim*
			      perm->outvec_leng);
       offset1 = 0;
       offset2 = 0;
       for (j = 0; j < ag->nullspace_dim; j++) {
	 ML_Operator_Apply(perm, perm->invec_leng, 
			   &((ag->nullspace_vect)[offset1]), 
			   perm->outvec_leng, &(new_null[offset2]));

	 offset1 += perm->invec_leng;
	 offset2 += perm->outvec_leng;
       }
       ML_Aggregate_Set_NullSpace(ag, ag->num_PDE_eqns, ag->nullspace_dim,
				  new_null,perm->outvec_leng);
       ML_free(new_null);

     }
     ML_Operator_Destroy(&perm);
     oldA = Amatrix;
     Amatrix= newA;
     newP = ML_Operator_Create(permt->comm);
     oldP = *Pmatrix;
     *Pmatrix = newP;
   }
#endif
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
/*ms*/
   default :
           if (mypid == 0) printf("ML_Aggregate_Coarsen : invalid scheme.\n");
           exit(1);
           break;
   } 
#ifdef ML_repartition

   /* restore Amatrix and delete repartitioned one */

   if (status == 0) {
     Amatrix = oldA;
     ML_Operator_Destroy(&newA);

     /* do a mat-mat mult to get the appropriate P for the */
     /* unpartitioned matrix.                              */

     ML_2matmult(permt, *Pmatrix, oldP, ML_CSR_MATRIX);
     ML_Operator_Destroy(Pmatrix);
     ML_Operator_Destroy(&permt);
     *Pmatrix = oldP;
   }
#endif


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
   if (comm->ML_mypid == 0)
      printf("Aggregation time \t= %e\n",t0);
#endif

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
   if (nbytes > 0) ML_memory_alloc( (void **) &request, nbytes, "AE1" );
   else            request = NULL;
   total_send_leng = 0;
   for ( i = 0; i < N_send_neighbors; i++ ) 
      total_send_leng += aggr_comm->send_leng[i];
   nbytes = total_send_leng * sizeof(double);
   if (nbytes > 0) ML_memory_alloc( (void **) &send_buf, nbytes, "AE2" );
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
      comm->USR_irecvbytes((void *) &(vec_data[offset]), length, &fromproc,
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
      comm->USR_sendbytes((void *) &(send_buf[offset]), length, toproc, 
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
      comm->USR_cheapwaitbytes((void *) &(vec_data[offset]), length, &fromproc,
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
   if ( nbytes > 0 ) ML_memory_alloc( (void **) &Request, nbytes, "AGZ" );
   else              Request = NULL;
   offset = 0;
   msgtype = msgid;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      fromproc = neighbors[i];
      length = recv_leng[i] * typeleng;
      comm->USR_irecvbytes(&recvbuf[offset*typeleng],length,&fromproc,
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
      comm->USR_sendbytes((void*) &sendbuf[offset*typeleng], length,
                          neighbors[i], msgtype, comm->USR_comm );
      offset += send_leng[i];
   }
   offset = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      fromproc = neighbors[i];
      length = recv_leng[i] * typeleng;
      msgtype = msgid;
      comm->USR_cheapwaitbytes(&recvbuf[offset*typeleng], length, &fromproc,
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
/* its role is to releave stress from aggregation routines whose quality     */
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
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ML_repartition_matrix(ML_Operator *mat, ML_Operator **new_mat,
			  ML_Operator **permutation, ML_Operator **permt,
			  int num_PDE_eqns)
			  
{
  int    proc_index, Nsnd, Nrcv, *neighbors, flag, oldj;
  int    *NeighborList, *Nsnds, *Nrcvs, **SndIndices, **RcvIndices, count;
  int    *the_list, the_length, offset, Nnonzero, Nglobal, oldNnonzero;
  int    *remote_offsets, Nneighbors, i, j, Nrows, Nghost;
  int    *permute_array, *itemp, *columns, *rowptr;
  double *dvec, *d2vec, *values;
  ML_Operator *perm_mat, *permt_mat, *permuted_Amat;
  ML_CommInfoOP *mat_comm = NULL;
  struct ML_CSR_MSRdata *temp;

  *new_mat = NULL;
  *permutation = NULL;
  Nglobal = mat->invec_leng;
  ML_gsum_scalar_int(&Nglobal, &i, mat->comm);

  if (Nglobal/num_PDE_eqns >= 12*mat->comm->ML_nprocs) return 1;

  /* Choose a random subset of global unknowns. These will become root nodes */
  /* for a very simple aggregation scheme. These aggregates/subdomains are   */ 
  /* later assigned to processors defining the new partitioning.             */

  ML_random_global_subset(mat, 20.,&the_list, &the_length, num_PDE_eqns);
  for (i = 0; i < the_length; i++) the_list[i] *= num_PDE_eqns;
#ifdef DEBUG
  if (mat->comm->ML_mypid == 0) {
    for (i = 0; i < the_length; i++) 
      printf("%d: ML_reparition_matrix: root_node(%d) = %d\n",mat->comm->ML_mypid,i,the_list[i]);
  }
  fflush(stdout);
#endif

  if (the_length >= mat->comm->ML_nprocs) {
    ML_free(the_list);
    return 1;
  }
  printf("%d: continue\n",mat->comm->ML_mypid); fflush(stdout);

  /* Make up a global numbering for the matrix unknowns. The jth  */
  /* unknown on processor k has global id = remote_offsets[k] + j */

  remote_offsets = (int *) ML_allocate(sizeof(int)*(mat->comm->ML_nprocs + 1));
  itemp         = (int *) ML_allocate(sizeof(int)*(Nglobal +
						   mat->comm->ML_nprocs + 1));
  for (i = 0; i < mat->comm->ML_nprocs; i++) remote_offsets[i] = 0;
  remote_offsets[mat->comm->ML_mypid] = mat->invec_leng;
  ML_gsum_vec_int(&remote_offsets,&itemp,mat->comm->ML_nprocs,mat->comm);
  j = 0; oldj = 0;
  for (i = 1; i < mat->comm->ML_nprocs; i++) {
    oldj = j;
    j += remote_offsets[i-1];
    remote_offsets[i-1] = oldj;
  }
  remote_offsets[ mat->comm->ML_nprocs-1] = j;
  remote_offsets[mat->comm->ML_nprocs] = Nglobal;
  offset = remote_offsets[mat->comm->ML_mypid];
#ifdef DEBUG
  if (mat->comm->ML_mypid == 0) 
    for (i = 0; i <= mat->comm->ML_nprocs; i++)
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
    ML_gsum_scalar_int(&Nnonzero,&i,mat->comm);
  }
  ML_free(d2vec);

  /* Assign any singletons to the a new domain */

  for (i = 0; i < mat->outvec_leng; i++) {
    if ( dvec[i] == 0.) dvec[i] = (double) (the_length + 1);
  }

#ifdef DEBUG
  for (i = 0; i < mat->outvec_leng; i++) 
    printf("%d: ML_repartition_matrix: domain_assignment(%d) = %e\n",mat->comm->ML_mypid,i,dvec[i]);
  fflush(stdout);
#endif

  /* Build a global array containing all the subdomain assignments */

  permute_array = (int *) ML_allocate(sizeof(int)*(Nglobal+1));
  for (i = 0; i <= Nglobal; i++) permute_array[i] = 0;
  for (i = 0; i < mat->invec_leng; i++) 
    permute_array[i+offset] = (int) dvec[i];
  ML_gsum_vec_int(&permute_array,&itemp,Nglobal,mat->comm);
  ML_free(itemp);
  ML_free(dvec);

  /* Count the number of rows assigned to this processor */
  Nrows = 0;
  for (i = 0; i < Nglobal; i++) {
    if (permute_array[i] - 1  == mat->comm->ML_mypid) Nrows++;
  }

  /* Now build a CSR matrix to represent the permuation */
  /* along with the communication information.          */

  columns = (int *) ML_allocate(sizeof(int)*Nrows);
  values  = (double *) ML_allocate(sizeof(double)*Nrows);
  rowptr = (int *) ML_allocate(sizeof(int)*(Nrows+1));
  neighbors       = (int *) ML_allocate(sizeof(int)*(Nglobal+1));
  for (i = 0; i < Nglobal; i++) neighbors[i] = -1;
  j = 0;	Nghost = 0; proc_index = 0; Nsnd = 0; Nrcv = 0;
  flag = 0;


  /* First determine with which processors do we need to */
  /* receive information to perform the permutation.     */

  NeighborList =  (int *) ML_allocate(sizeof(int)*
				      ML_min(Nglobal+1, mat->comm->ML_nprocs));
  Nneighbors = 0;
  proc_index = 0;
  for (i = 0; i <= Nglobal; i++) {
    if (i != Nglobal) {
      while (i == remote_offsets[proc_index+1]) {
	proc_index++;
      }
    }
    if ( (permute_array[i] - 1 == mat->comm->ML_mypid) &&
	 ( (i < offset) || (i >= offset + mat->invec_leng) )) {
      if (neighbors[proc_index] == -1) {
	neighbors[proc_index] = Nneighbors;
	NeighborList[Nneighbors++] = proc_index;
      }
    }
  }
	    
  /* Where (which processors) do we need to send information */

  for (i = offset; i < offset + mat->invec_leng; i++) {
    if (permute_array[i] - 1 != mat->comm->ML_mypid) {
      if (neighbors[permute_array[i]-1] == -1) {
	neighbors[permute_array[i]-1] = Nneighbors;
	NeighborList[Nneighbors++] = permute_array[i]-1;
      }
    }
  }

  /* Count how much information is sent and received from */
  /* each processor in the NeighborList.                  */

  Nsnds = (int *) ML_allocate(sizeof(int)*(Nneighbors+1));
  Nrcvs = (int *) ML_allocate(sizeof(int)*(Nneighbors+1));
  for (i = 0; i < Nneighbors; i++) Nsnds[i] = 0;
  for (i = 0; i < Nneighbors; i++) Nrcvs[i] = 0;

  for (i = offset; i < offset + mat->invec_leng; i++) {
    if (permute_array[i] - 1 != mat->comm->ML_mypid)
      Nsnds[neighbors[permute_array[i]-1]]++;
  }

  proc_index = 0;
  for (i = 0; i <= Nglobal; i++) {
    if (i != Nglobal) {
      while (i == remote_offsets[proc_index+1]) proc_index++;
    }
    if ( (permute_array[i] - 1 == mat->comm->ML_mypid) &&
	 ( (i < offset) || (i >= offset + mat->invec_leng) )) {
      Nrcvs[neighbors[proc_index]]++;
    }
  }

  /* Now fill in the indices which must be sent and received */
  SndIndices = (int **) ML_allocate(sizeof(int)*(Nneighbors+1));
  RcvIndices = (int **) ML_allocate(sizeof(int)*(Nneighbors+1));
  for (i = 0; i < Nneighbors; i++) {
    SndIndices[i] = (int *) ML_allocate(sizeof(int)*(Nsnds[i]+1));
    RcvIndices[i] = (int *) ML_allocate(sizeof(int)*(Nrcvs[i]+1));
    Nsnds[i] = 0;
    Nrcvs[i] = 0;
  }

  for (i = offset; i < offset + mat->invec_leng; i++) {
    if (permute_array[i] - 1 != mat->comm->ML_mypid) {
      j = neighbors[permute_array[i]-1];
      SndIndices[j][Nsnds[j]] = i - offset;
      Nsnds[j]++;
    }
  }

  count = mat->invec_leng;
  proc_index = 0;
  for (i = 0; i <= Nglobal; i++) {
    if (i != Nglobal) {
      while (i == remote_offsets[proc_index+1]) proc_index++;
    }
    if ( (permute_array[i] - 1 == mat->comm->ML_mypid) &&
	 ( (i < offset) || (i >= offset + mat->invec_leng) )) {
      j = neighbors[proc_index];
      RcvIndices[j][Nrcvs[j]] = count++;
      Nrcvs[j]++;
    }
  }
  ML_free(remote_offsets);
  ML_free(neighbors);


  ML_CommInfoOP_Set_neighbors(&mat_comm, Nneighbors, 
			      NeighborList,ML_OVERWRITE, NULL, 0);


  for (i = 0; i < Nneighbors; i++) {
    ML_CommInfoOP_Set_exch_info(mat_comm, NeighborList[i],
				Nrcvs[i], RcvIndices[i], Nsnds[i], SndIndices[i]);
    ML_free(SndIndices[i]);
    ML_free(RcvIndices[i]);
  }
  ML_free(NeighborList);
  ML_free(Nsnds);
  ML_free(Nrcvs);
  ML_free(SndIndices);
  ML_free(RcvIndices);


  /* Now put in the CSR data */

  j = 0;	Nghost = 0; proc_index = 0; Nsnd = 0; Nrcv = 0;
  for (i = 0; i < Nglobal; i++) {
    if (permute_array[i] - 1 == mat->comm->ML_mypid) {
      rowptr[j] = j;
      values[j] = 1.;
      if ( (i >= offset) && 
	   (i < offset + mat->invec_leng) ) {
	/* local row */
	columns[j] = i - offset;
      }
      else {
	/* remote row */
	columns[j] = mat->invec_leng + Nghost;
	Nghost++;
      }
      j++;
    }
  }
  ML_free(permute_array);
  rowptr[Nrows] = Nrows;
  perm_mat = ML_Operator_Create(mat->comm);
  temp = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  temp->columns = columns;
  temp->values = values;
  temp->rowptr = rowptr;

  ML_Operator_Set_ApplyFuncData(perm_mat,mat->invec_leng,
				Nrows,ML_EMPTY,temp,mat->invec_leng,NULL,0);
  ML_Operator_Set_Getrow(perm_mat, ML_EXTERNAL, Nrows, CSR_getrows);
  ML_Operator_Set_ApplyFunc( perm_mat, ML_INTERNAL, CSR_matvec);
  perm_mat->getrow->pre_comm = mat_comm;
  perm_mat->data_destroy = ML_CSR_MSRdata_Destroy;
  permt_mat = ML_Operator_Create(mat->comm);

  ML_Operator_Transpose_byrow(perm_mat,permt_mat);

  permuted_Amat = ML_Operator_Create(mat->comm);
  ML_rap(perm_mat,mat,permt_mat,permuted_Amat, ML_CSR_MATRIX);
  *new_mat = permuted_Amat;
  *permutation = perm_mat;
  *permt       = permt_mat;

  return 0;
}

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




