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

/* ************************************************************************* */
/* Constructor                                                               */
/* ------------------------------------------------------------------------- */
#ifdef AZTEC
extern int ML_Aggregate_AztecRead(ML_Aggregate *ag);
#endif

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
   (*ag)->threshold                  = 0.08;
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
   ML_memory_alloc((void*) &(ag->aggr_info), level*sizeof(int*),"AGu");
   for ( i = 0; i < level; i++ ) ag->aggr_info[i] = NULL;
   ML_memory_alloc((void*) &(ag->aggr_count), level*sizeof(int),"AGx");
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
   if ( level < 0 || level >= ag->max_levels ) 
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
   int i=1, ndofs, Ncoarse, coarsen_scheme;
   int mypid, nprocs;
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
      printf("ML_Aggregate_Coarsen begins \n");

/* #### moved this somewhere else ?? */
   Amatrix->num_PDEs = ag->num_PDE_eqns;
   Amatrix->num_rigid = ag->nullspace_dim;

   /* check to see which aggregation algorithm to use */

   ndofs = Amatrix->outvec_leng;
   if ( ndofs < 2 ) ndofs = 0; else ndofs = 1;
   ML_gsum_vec_int(&ndofs, &i, 1, comm);
   /*
   gmin = ML_gmax_double((double) (-1.0 * Amatrix->outvec_leng) , comm);
   if (comm->ML_mypid == 0)
      printf("Smallest Amatrix->outvec_leng = %d\n", (int) (-1 * gmin));
   */
   nprocs = comm->ML_nprocs;
   coarsen_scheme = ag->coarsen_scheme;
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
      else
      {
         (*Pmatrix) = NULL;
         return 0;
      }
   }
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
   if (comm->ML_mypid == 0)
      printf(" Aggregation time \t= %e\n",t0);
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
   switch (ag->coarsen_scheme)
   {
      case ML_AGGR_UNCOUPLED :
           printf("ML_Aggregate : coarsen scheme     = UNCOUPLED\n");
           break;
      case ML_AGGR_COUPLED :
           printf("ML_Aggregate : coarsen scheme     = COUPLED\n");
           break;
   }
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
      comm->USR_waitbytes((void *) &(vec_data[offset]), length, &fromproc,
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
                        &msgtype, comm->USR_comm, (void *) &Request[i] );
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
      comm->USR_waitbytes(&recvbuf[offset*typeleng], length, &fromproc,
                          &msgtype, comm->USR_comm, (void *) &Request[i] );
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


