/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* ******************************************************************** */
/* Functions to create tentative prolongators                           */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : December, 1999                                       */
/* ******************************************************************** */
/* ******************************************************************** */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"

#ifdef ML_AGGR_READINFO
#include "az_aztec.h"
#endif

#define dabs(x) (((x) > 0) ? x : (-(x)))

/* ******************************************************************** */
/* variables used for parallel debugging  (Ray)                         */
/* -------------------------------------------------------------------- */

#ifdef ML_AGGR_PARTEST
extern int **global_mapping = NULL, global_nrows, global_ncoarse;
#endif

/* ******************************************************************** */
/* internal function defined later on in this file                      */
/* -------------------------------------------------------------------- */

extern void ML_CSR_MSR_ML_memorydata_Destroy(void *data);

/* ******************************************************************** */
/* external functions called from this file                             */
/* -------------------------------------------------------------------- */

extern int ML_randomize(int nlist, int *list);
extern int ML_Aggregate_CoarsenUncoupledVBlock(ML_Aggregate*,ML_Operator*, 
                                               ML_Operator **, ML_Comm *);
extern int ML_Aggregate_CoarsenCoupledVBlock(ML_Aggregate*,ML_Operator*, 
                                             ML_Operator **, ML_Comm *);

/* ******************************************************************** */
/* local defines                                                        */
/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
/* Constructor                                                          */
/* -------------------------------------------------------------------- */
int phase2_off = 0;

int ML_Aggregate_Create( ML_Aggregate **ag )
{
#ifdef ML_AGGR_READINFO
   FILE *fp;
   int  proc_config[AZ_PROC_SIZE];
#endif
	 
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
   (*ag)->max_levels                 = 0;
   (*ag)->max_coarse_size            = 100;
   (*ag)->begin_level                = 0;
   (*ag)->cur_level                  = 0;
   (*ag)->aggr_info                  = NULL;
   (*ag)->aggr_count                 = NULL;
   (*ag)->drop_tol_for_smoothing     = 0.0;
   (*ag)->fine_complexity            = 0.0;
   (*ag)->nvblocks                   = 0;
   (*ag)->vblock_info                = NULL;
   (*ag)->operator_complexity        = 0.0;

#ifdef ML_AGGR_READINFO
#ifdef ML_MPI
   AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
   AZ_set_proc_config(proc_config, AZ_NOT_MPI );
#endif

   if (proc_config[AZ_node] == 0) {
      fp = fopen("PaRams","r");
      if (fp == NULL) { printf("woops no PaRams file\n"); exit(1);}
      fscanf(fp,"%d", &((*ag)->ordering) );
      fscanf(fp,"%d", &((*ag)->min_nodes_per_aggregate) );
      fscanf(fp,"%d", &((*ag)->max_neigh_already_selected) );
      fscanf(fp,"%d", &((*ag)->attach_scheme) );
      fscanf(fp,"%d", &((*ag)->max_levels) );
      fscanf(fp,"%d", &((*ag)->coarsen_scheme) );
      fscanf(fp,"%lf", &((*ag)->threshold) );
      fscanf(fp,"%lf", &((*ag)->smoothP_damping_factor) );
      fscanf(fp,"%lf", &((*ag)->drop_tol_for_smoothing) );
      fscanf(fp,"%d", &phase2_off );
      fclose(fp);
    }
    AZ_broadcast((char*)&((*ag)->ordering),sizeof(int),proc_config,AZ_PACK);
    AZ_broadcast((char*)&((*ag)->min_nodes_per_aggregate),sizeof(int), 
                  proc_config, AZ_PACK);
    AZ_broadcast((char*)&((*ag)->max_neigh_already_selected),sizeof(int), 
                  proc_config, AZ_PACK);
    AZ_broadcast((char*)&((*ag)->attach_scheme),sizeof(int),proc_config,
                  AZ_PACK);
    AZ_broadcast((char*)&((*ag)->max_levels),sizeof(int),proc_config,AZ_PACK);
    AZ_broadcast((char*)&((*ag)->coarsen_scheme),sizeof(int),proc_config,
                  AZ_PACK);
    AZ_broadcast((char*)&((*ag)->threshold),sizeof(double),proc_config,AZ_PACK);
    AZ_broadcast((char*)&((*ag)->smoothP_damping_factor), sizeof(double), 
                  proc_config, AZ_PACK);
    AZ_broadcast((char*)&((*ag)->drop_tol_for_smoothing), sizeof(double), 
                  proc_config, AZ_PACK);
    AZ_broadcast((char *) &(phase2_off), sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char*)NULL         , 0          , proc_config, AZ_SEND);
#endif

   return 0;
}

/* ******************************************************************** */
/* destructor                                                           */
/* -------------------------------------------------------------------- */

int ML_Aggregate_Destroy( ML_Aggregate **ag )
{
   int i;

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
   ML_memory_free( (void **) ag );
   (*ag) = NULL;
   return 0;
}

/* ******************************************************************** */
/* set minimum number of nodes per aggregate for phase 2                */
/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
/* set output level                                                     */
/* -------------------------------------------------------------------- */

int ML_Aggregate_Set_OutputLevel( ML_Aggregate *ag, int level )
{
   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      printf("ML_Aggregate_Set_OutputLevel : wrong object. \n");
      exit(-1);
   }
   ag->print_flag = level;
   return 0;
}

/* ******************************************************************** */
/* set random or natural ordering for traversing the matrix graph       */
/* -------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
/* select scheme to put un-aggregated nodes into existing aggregates    */
/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
/* set maximum coarsest grid size                                       */
/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
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

/* -------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
/* set/reset aggregation threshold                                      */
/* -------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
/* set damping factor for the smoothed prolongator                      */
/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
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

/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
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

/* ******************************************************************** */
/* set max number of levels for the smoothed prolongator                */
/* -------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
/* get the aggregation information for a given level                    */
/* -------------------------------------------------------------------- */

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

/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
/* set nullspace                                                        */
/* null_vect can be NULL, in which case the default nullspace           */
/* (ones for each PDE) will be used.  leng is the dimension of A at the */
/* coarsest level                                                       */
/* -------------------------------------------------------------------- */

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

   if (ag->nullspace_vect != NULL)
      ML_memory_free((void **)&(ag->nullspace_vect));

   /* if the user-supplied nullspace vector isn't null, allocate space */
   /* and load it */

#ifdef ML_AGGR_NSOUTPUT
   sprintf( fname, "null.%d", global_comm->ML_mypid);
   fp = fopen(fname, "w");
#endif

   if (null_vect != NULL) 
   {
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

/* ******************************************************************** */
/* ******************************************************************** */
/* Coarsening routine                                                   */
/* -------------------------------------------------------------------- */

int ML_Aggregate_Coarsen( ML_Aggregate *ag, ML_Operator *Amatrix, 
                          ML_Operator **Pmatrix, ML_Comm *comm)
{
   int i, ndofs, Ncoarse, coarsen_scheme;
   int mypid;
#ifdef ML_TIMING
   double t0;
#endif
#ifdef ML_TIMING
   t0 = GetClock();
#endif

   mypid = comm->ML_mypid;

   if ( ag->ML_id != ML_ID_AGGRE ) 
   {
      printf("ML_Aggregate_Coarsen : wrong object. \n");
      exit(-1);
   }

   if ( mypid == 0 && ag->print_flag ) 
      printf("ML_Aggregate_Coarsen begins ...\n");

/* #### moved this somewhere else ?? */
   Amatrix->num_PDEs = ag->num_PDE_eqns;
   Amatrix->num_rigid = ag->nullspace_dim;

   ndofs = Amatrix->outvec_leng;
   if ( ndofs < 2 ) ndofs = 0; else ndofs = 1;
   i = 1;
   ML_gsum_vec_int(&ndofs, &i, 1, comm);
   coarsen_scheme = ag->coarsen_scheme;

   if ( coarsen_scheme == ML_AGGR_UNCOUPLED )
   {
      /*-------------------------------------------------------------- */
      /* This portion of the code is for testing variable dof per node */
      /*                                                               */
      /*if ((ag->begin_level-ag->cur_level)==0)                        */
      /*{                                                              */
      /*   ag->nvblocks = Amatrix->outvec_leng/4*3/2;                  */
      /*   ag->vblock_info = (int *) malloc(ag->nvblocks*sizeof(int)); */
      /*   for (i=0; i < ag->nvblocks; i+=3 )                          */
      /*   {                                                           */
      /*      ag->vblock_info[i] = 4;                                  */
      /*      ag->vblock_info[i+1] = 2;                                */
      /*      ag->vblock_info[i+2] = 2;                                */
      /*   }                                                           */
      /*}                                                              */
      /*-------------------------------------------------------------- */

      ag->nvblocks = Amatrix->outvec_leng / ag->num_PDE_eqns;
      ag->vblock_info = (int *) malloc(ag->nvblocks*sizeof(int));
      for (i=0; i < ag->nvblocks; i++ ) ag->vblock_info[i] = ag->num_PDE_eqns;

      /* Cant find VBlock so I am taking it out */

      Ncoarse = ML_Aggregate_CoarsenUncoupled(ag,Amatrix,
                                                    Pmatrix,comm);
      ag->nvblocks = 0;
      free(ag->vblock_info);

      /* rst: I can't find Vblock routines any more so I took out the */
      /* invokation. */

      /*-------------------------------------------------------------- 
      if ((ag->max_levels-ag->cur_level-1)==0 && ag->nvblocks>0)
         Ncoarse=ML_Aggregate_CoarsenUncoupled(ag,Amatrix,Pmatrix,comm);
      else
         Ncoarse = ML_Aggregate_CoarsenUncoupled(ag,Amatrix,Pmatrix,comm);
       *-------------------------------------------------------------- */
   } 
   else if ( coarsen_scheme == ML_AGGR_COUPLED )
   {
      ag->nvblocks = Amatrix->outvec_leng / ag->num_PDE_eqns;
      ag->vblock_info = (int *) malloc(ag->nvblocks*sizeof(int));
      for (i=0; i < ag->nvblocks; i++ ) ag->vblock_info[i] = ag->num_PDE_eqns;
/*
      if ((ag->begin_level-ag->cur_level)==0)
      {
         ag->nvblocks = Amatrix->outvec_leng/2;
         if ( ag->nvblocks*2 != Amatrix->outvec_leng)
            printf("WARNING : blocks cannot be formed.\n");
         ag->vblock_info = (int *) malloc(ag->nvblocks*sizeof(int));
         for (i=0; i < ag->nvblocks; i++ )
         {
            ag->vblock_info[i] = 2;
         }
      }
      else
*/
      {
         ag->nvblocks = Amatrix->outvec_leng;
         ag->vblock_info = (int *) malloc(ag->nvblocks*sizeof(int));
         for (i=0; i < ag->nvblocks; i++ ) ag->vblock_info[i] = 1;
      }
/*
Ncoarse = ML_Aggregate_CoarsenCoupled(ag,Amatrix,Pmatrix,comm);
*/
      /* rst: I can't find Vblock routines any more so I took out the */
      /* invokation. */

      Ncoarse = ML_Aggregate_CoarsenCoupled(ag,Amatrix,Pmatrix,comm);
      ag->nvblocks = 0;
      free(ag->vblock_info);

      /*-------------------------------------------------------------- */
      /*ag->nvblocks = Amatrix->outvec_leng;                           */
      /*ag->vblock_info = (int *) malloc(ag->nvblocks*sizeof(int));    */
      /*for (i=0; i < ag->nvblocks; i++ ) ag->vblock_info[i] = 1;      */
      /*                                                               */
      /*Ncoarse=ML_Aggregate_CoarsenCoupledVBlock(ag,Amatrix,Pmatrix,comm);*/
      /*                                                               */
      /*ag->nvblocks = 0;                                              */
      /*free(ag->vblock_info);                                         */
      /*-------------------------------------------------------------- */
   }
   else if ( coarsen_scheme == ML_AGGR_MIS )
   {
      Ncoarse = ML_Aggregate_CoarsenMIS(ag,Amatrix,Pmatrix,comm);
   }
   else if ( coarsen_scheme == ML_AGGR_DD )
   {
      Ncoarse = ML_Aggregate_CoarsenDomainDecomp(ag,Amatrix,Pmatrix,comm);
   }
   else 
   {
      if ( mypid == 0 ) printf("ML_Aggregate_Coarsen : invalid scheme.\n");
      exit(1);
   } 

#ifdef ML_AGGR_DEBUG
   i = 0;
   i = ML_gmax_int(i, comm);
   if ( mypid == 0 && ag->print_flag ) 
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

/* ******************************************************************** */
/* ******************************************************************** */
/* construct the tentative prolongator (local)                          */
/*  phase 1 : relax on the new seed point as Vanek                      */
/*  phase 2 : assign the rest of the nodes to one of the existing       */
/*            aggregate (attach_scheme), if possible.                   */
/*  phase 3 : see if the un-aggregated nodes have enough neighbors      */
/*            (min_nodes_per_aggregate) to form its own aggregate       */
/* -------------------------------------------------------------------- */

int ML_Aggregate_CoarsenUncoupled(ML_Aggregate *ml_ag,ML_Operator *Amatrix, 
                                  ML_Operator **Pmatrix, ML_Comm *comm)
{
   int     i, j, k, m, inode = 0, jnode, nbytes, length, *mat_indx=NULL, Nrows;
   int     select_flag, aggr_count, index, mypid, Ncoarse, inode2;
   int     *aggr_index = NULL, search_flag, *itmp_array = NULL, count;
   int     mincount, *int_buf = NULL, *aggr_stat = NULL, diff_level;
   int     *new_ia = NULL, *new_ja = NULL, maxcount, *randomVector = NULL;
   int     aggr_cnt_leng, *aggr_cnt_array=NULL, maxnnz_per_row=500;
   int     min_nodes_per_aggregate, max_neigh_selected, attach_scheme;
   int     ordering, lwork, *agg_sizes=NULL, num_PDE_eqns;
   int     *agg_sizes_cum = NULL, nz_cnt, *amal_mat_indx=NULL, amal_count;
   int     kk, **rows_in_aggs = NULL, max_agg_size, info, level, print_flag;
   int     ii,jj, block_col, row, nullspace_dim;
   char    *col_entered;
   double  dcompare1, dcompare2, *new_val = NULL, *diagonal = NULL, epsilon;
   double  *dble_buf = NULL, *tmp_vect = NULL, epsmax;
   double  *nullspace_vect=NULL, *new_null=NULL, *work=NULL, *qr_tmp=NULL;
   struct  ML_CSR_MSRdata *csr_data;
   ML_Node           *node_head = NULL, *node_tail = NULL, *new_node;
   ML_SuperNode      *aggr_head, *aggr_curr = NULL, *supernode;
   ML_Aggregate_Comm *aggr_comm;
   ML_GetrowFunc     *getrow_obj;
   int               (*getrowfunc)(void *,int,int*,int,int*,double*,int*);

#ifdef ML_AGGR_TEST
   int     nprocs, *itmp_array2, Nrows_offset, naggr_offset;
   char    zfn[80];
   FILE    *zfp;
#endif

#ifdef ML_AGGR_PARTEST
   int   proc, globrow, cc;
   int   local_to_global_row(int i);
   FILE  *pfp, *nfp, *afp;
   char  pfn[80], nfn[80], afn[80];
#endif

#ifdef MATLAB2
   FILE  *afp;
   char  afn[80];
#endif
double thesign, largest;

   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   mypid                   = comm->ML_mypid;
#ifdef ML_AGGR_TEST
   nprocs                  = comm->ML_nprocs;
#endif
   min_nodes_per_aggregate = ml_ag->min_nodes_per_aggregate;
   max_neigh_selected      = ml_ag->max_neigh_already_selected;
   epsilon                 = ml_ag->threshold;
   ordering                = ml_ag->ordering;
   attach_scheme           = ml_ag->attach_scheme;
   num_PDE_eqns            = ml_ag->num_PDE_eqns;
   nullspace_dim           = ml_ag->nullspace_dim;
   nullspace_vect          = ml_ag->nullspace_vect;
   Nrows                   = Amatrix->outvec_leng;

   /* ============================================================= */
   /* check the system size versus null dimension size              */
   /* ============================================================= */

   /*if ( num_PDE_eqns != nullspace_dim && nullspace_vect == NULL )
    *{
    *   nullspace_dim = num_PDE_eqns;
    *   ml_ag->nullspace_dim = nullspace_dim;
    *} */ 
   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level < 0 ) diff_level = - diff_level;
   if ( diff_level != 0 ) num_PDE_eqns = nullspace_dim; /* ## 6/14/00 */

   /* ============================================================= */
   /* set up the threshold for weight-based coarsening              */
   /* ============================================================= */

   epsmax = 0.25;
   while ( diff_level > 0 ) {epsmax *= 0.5; epsilon *= 0.5; diff_level--;}
   /*if (epsilon < epsmax && epsilon != 0.0) 
    *   epsilon += (epsmax - epsilon) / 100 * nprocs;
    *if (epsilon > epsmax) epsilon = epsmax; */

#ifdef ML_AGGR_DEBUG
   if ( mypid == 0 )
   {
      printf("ML_Aggregate_CoarsenUncoupled : current level = %d\n", 
                                              ml_ag->cur_level);
      printf("ML_Aggregate_CoarsenUncoupled : current eps = %e\n",epsilon);
      printf("ML_Aggregate_CoarsenUncoupled : block size  = %d\n",num_PDE_eqns);
   }
#endif
   epsilon = epsilon * epsilon;

   /* ============================================================= */
   /* Construct an initial status array to store whether the nodes  */
   /* have been aggregated. aggr_index stores the aggregate number  */
   /* where this node has been aggregated into.                     */
   /* ============================================================= */

   nbytes = Nrows / num_PDE_eqns * sizeof( int );
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &aggr_index, nbytes, "AGG");
      ML_memory_alloc((void**) &aggr_stat,  nbytes, "AGH");
   } else aggr_index = aggr_stat = NULL;
   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) aggr_stat[i] = ML_AGGR_READY;
   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) aggr_index[i] = -1;

   /* ============================================================= */
   /* generate MSR matrix from incoming A matrix                    */
   /* ============================================================= */

   /* ------------------------------------------------------------- */
   /* first find out whether the getrow function is available       */
   /* ------------------------------------------------------------- */

   getrow_obj = Amatrix->getrow;
   if ( getrow_obj->ML_id == ML_EXTERNAL) 
   {
      getrowfunc = getrow_obj->external;
   } 
   else 
   {
      getrowfunc = getrow_obj->internal;
   }
   if ( getrowfunc == NULL ) 
   {
      printf("ML_Aggregate_Coarsen Error : no getrow function.\n");
      exit(-1);
   }

   /* ------------------------------------------------------------- */
   /* allocate initial temporary storage space for getrow           */
   /* also allocate space for storing the diagonal (if epsilon>0)   */
   /* ------------------------------------------------------------- */

   nbytes = maxnnz_per_row * sizeof( int );
   ML_memory_alloc((void**) &int_buf, nbytes, "AGA");
   nbytes = maxnnz_per_row * sizeof( double );
   ML_memory_alloc((void**) &dble_buf, nbytes, "AGB");
   if ( Nrows > 0 ) 
   {
      nbytes = Nrows * sizeof( double );
      ML_memory_alloc((void**) &diagonal, nbytes, "AGC");
   } else diagonal = NULL;
   if ( Nrows / num_PDE_eqns * num_PDE_eqns != Nrows )
   {
      printf("ML_Aggregate_Coarsen Error : incomplete blocks %d %d.\n",
              Nrows, num_PDE_eqns);
      exit(-1);
   }

   /* ------------------------------------------------------------- */
   /* find out about how much memory to allocate for the matrix     */
   /* ------------------------------------------------------------- */

   count = 0;
   for ( i = 0; i < Nrows; i++ ) 
   {
      diagonal[i] = 0.0;
      while (getrowfunc(Amatrix->data, 1, &i, maxnnz_per_row, int_buf, 
                          dble_buf, &m) == 0 ) 
      {
         ML_memory_free((void**) &int_buf);
         ML_memory_free((void**) &dble_buf);
         maxnnz_per_row = maxnnz_per_row * 2 + 1; 
         nbytes = maxnnz_per_row * sizeof( int );
         ML_memory_alloc((void**) &int_buf, nbytes, "AGD");
         nbytes = maxnnz_per_row * sizeof( double );
         ML_memory_alloc((void**) &dble_buf,  nbytes, "AGE");
      }
      for ( j = 0; j < m; j++ ) 
      {
         if ( int_buf[j] == i ) diagonal[i] = dble_buf[j];
         /*if (dabs(dble_buf[j]) > diagonal[i] ) diagonal[i] = dabs(dble_buf[j]);*/
      }
      count += m;
      if ( diagonal[i] == 0.0 ) count++;
   }
   print_flag = 1;
   for ( i = 0; i < Nrows; i++ ) 
   {
      if ( diagonal[i] == 0.0 )
      {
         if ( print_flag == 1 ) 
         {
            print_flag = 0;
            printf("Aggregation Coarsening warning : some diag = 0.0\n");
         }
         diagonal[i] = 1.0;
      }
   }
   if ( epsilon == 0.0 && diagonal != NULL ) 
   {
      ML_memory_free((void**) &diagonal);
      diagonal = NULL;
   }

   /* ------------------------------------------------------------- */
   /* allocate memory for the entire matrix (only the column indices*/
   /* are needed since the matrix will be pruned here               */
   /* ------------------------------------------------------------- */

   nbytes = (count + 1) * sizeof( int );
   ML_memory_alloc((void**) &mat_indx, nbytes, "AGF");
   k = ML_Comm_GsumInt( comm, Nrows);
   m = ML_Comm_GsumInt( comm, count);

#ifdef ML_AGGR_DEBUG
   if ( mypid == 0 ) 
      printf("Aggregation(U) : Total nonzeros = %d (Nrows=%d)\n",m,k);
#endif   
   if ( ml_ag->operator_complexity == 0.0 )
   {
      ml_ag->fine_complexity = 1.0 * m;
      ml_ag->operator_complexity = 1.0 * m;
   }
   else
   {
      ml_ag->operator_complexity += 1.0 * m;
   }

   count = Nrows + 1;
   mat_indx[0] = count; 

   /* ------------------------------------------------------------- */
   /* extract the matrix using the getrow function                  */
   /* (pruning is done at this stage)                               */
   /* ------------------------------------------------------------- */

#ifdef ML_AGGR_PARTEST
   sprintf(afn,"AmatrixProc%dNc%d.dat",mypid,global_ncoarse);
   afp=fopen(afn,"w");
#endif
#ifdef MATLAB2
   sprintf(afn,"Amatrix2.dat");
   afp=fopen(afn,"w");
#endif

   for ( i = 0; i < Nrows; i++ ) 
   {
      getrowfunc(Amatrix->data,1,&i,maxnnz_per_row,int_buf,dble_buf, &m);
      if ( m > maxnnz_per_row ) printf("BIG WARNING\n");

      nz_cnt = 0;
      for (j = 0; j < m; j++) 
      {
#ifdef ML_AGGR_PARTEST
      /* now the diagonal entry locations (of the 4's) are okay, as
         are the elements in the same grid row (at positions +/- 2
         of the diag elts) but the vertical neighbors are just totally
         off with no obvious pattern.  How are the columns numbered?! */

         globrow=local_to_global_row(i);
         fprintf(afp,"%d %d %f\n", globrow, local_to_global_row(int_buf[j]),
                                dble_buf[j]);
#endif
#ifdef MATLAB2
         fprintf(afp,"%d %d %f\n", i, int_buf[j],dble_buf[j]);
#endif

         jnode = int_buf[j];
	 if ( jnode != i && jnode < Nrows && epsilon > 0.0 ) 
         {
            dcompare1 = dble_buf[j] * dble_buf[j];
            if ( dcompare1 > 0.0 )
            {
               nz_cnt++;
	       dcompare2 = diagonal[i] * diagonal[jnode];
	       if ( dcompare1 < 0 ) dcompare1 = - dcompare1;
	       if ( dcompare2 < 0 ) dcompare2 = - dcompare2;
	       if ( dcompare1 >= epsilon * dcompare2 ) 
	          mat_indx[count++] = int_buf[j];
            }
         } 
         else if ( jnode != i && jnode < Nrows && dble_buf[j] != 0.0)
         {
            mat_indx[count++] = int_buf[j];
            nz_cnt++;
         }
      }
      if ( nz_cnt == 0 ) {aggr_stat[i/num_PDE_eqns] = ML_AGGR_BDRY;}
      mat_indx[i+1] = count;
   }
   nz_cnt = count;

   /* ensure at least one aggregate per processor */

   count = 0;
   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) 
   {
      if ( aggr_stat[i] == ML_AGGR_READY ) count++;
   }
   if ( count == 0 )
   {
      for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) 
      {
         if ( aggr_stat[i] == ML_AGGR_BDRY ) 
         {
            aggr_stat[i] = ML_AGGR_READY;
            break;
         }
      }
   }

   ML_memory_free((void**) &int_buf);
   ML_memory_free((void**) &dble_buf);
   if ( diagonal != NULL ) ML_memory_free((void**) &diagonal);

#ifdef MATLAB
   fclose(afp);
   if ( mypid == 0 ) 
      printf("Aggregation (U) : Matrix read in\n");
#endif
#ifdef MATLAB2
   fclose(afp);
   if ( mypid == 0 ) 
      printf("Aggregation (U) : Matrix read in\n");
#endif

   /* ============================================================= */
   /* Construct the matrix that relates to the nodes by combining   */
   /* the rows of the matrix corresponding to different PDE's at    */
   /* the same node.                                                */
   /* ============================================================= */

   if (num_PDE_eqns == 1)  /* don't need to do anything */
   {
      amal_mat_indx = mat_indx;
   }
   else 
   {
      /* start by sorting the columns for each row of mat_indx to make it 
         easier to amalgamate the matrix, and initialize where we are in 
         each row */

#ifdef ML_AGGR_OBSOLETE
      ML_memory_alloc( (void**) &amal_locs, Nrows*sizeof(int), "MLs");

      for (i = 0; i < Nrows; i++) 
      {
         tmp = ML_sort(mat_indx[i+1]-mat_indx[i], mat_indx+mat_indx[i]);
         amal_locs[i] = mat_indx[i];
      }
#endif

      nbytes = (nz_cnt + 1) * sizeof( int ); /* probably excessive */
      if (nbytes > 0) ML_memory_alloc((void**) &amal_mat_indx,nbytes,"AGZ");

      amal_count = Nrows / num_PDE_eqns + 1;
      amal_mat_indx[0] = amal_count; 
      row = 0;
      col_entered = (char *) malloc(sizeof(char)*(1+ Nrows/num_PDE_eqns) );
      if (col_entered == NULL) {
         printf("Not enough space in ML_aggregate\n");
         exit(1);
      }
      for ( ii = 0; ii < Nrows/num_PDE_eqns; ii++) col_entered[ii] = 'F';

      for ( ii = 0; ii < Nrows/num_PDE_eqns; ii++) {
         col_entered[ii] = 'T';
         for ( jj = 0; jj < num_PDE_eqns; jj++) {
            for ( kk = mat_indx[row]; kk < mat_indx[row+1]; kk++) {
              block_col = mat_indx[kk]/num_PDE_eqns;
              if (col_entered[block_col] == 'F') {
                 amal_mat_indx[ amal_count++] = block_col;
                 col_entered[block_col] = 'T';
              }
            }
            row++;
         }
         amal_mat_indx[ii+1] = amal_count;
         col_entered[ii] = 'F';
         for ( jj = amal_mat_indx[ii]; jj < amal_mat_indx[ii+1]; jj++)
            col_entered[ amal_mat_indx[jj]] = 'F';
      }
      free(col_entered);

#ifdef ML_AGGR_OBSOLETE
      for ( i = 0; i < Nrows/num_PDE_eqns; i++) 
      {
         for (j = 0; j < Nrows/num_PDE_eqns; j++) 
         {
            if ( i !=j ) 
            { /* we assume the diagonal is non-zero, so we just have to */
              /*deal with off-diagonal elements */
               chng=0;
               for (kk=0; kk < num_PDE_eqns; kk++) 
               {
                  tmprow = i * num_PDE_eqns + kk; /* look at the row of A */
                  while((amal_locs[tmprow] < mat_indx[tmprow+1]) && 
                        (mat_indx[amal_locs[tmprow]] < (j+1)*num_PDE_eqns))
                  {
                     if (mat_indx[amal_locs[tmprow]] >= j*num_PDE_eqns) 
                        if (chng==0) 
                        {
                           amal_mat_indx[amal_count++]=j;
                           chng=1;
                        }
                     amal_locs[tmprow]++;
                  }
               }
               /*####### Tong (deleted) */
               /*if (chng==1) amal_mat_indx[i+1]=amal_count; */
            }
         }
         /*####### Tong (added) */
         amal_mat_indx[i+1] = amal_count;
      }
      ML_memory_free( (void**) &amal_locs);
#endif
   }

#ifdef ML_AGGR_DEBUG
   if (( mypid == 0 ) && (num_PDE_eqns > 1))
      printf("Aggregation (U) : Amalgamated matrix done \n");
#endif

   /* ============================================================= */
   /* Set up the data structures for aggregation                    */
   /* ============================================================= */

   aggr_count = 0;
   aggr_head = NULL;
   aggr_cnt_leng = Nrows / num_PDE_eqns / 5 + 2;
   nbytes = aggr_cnt_leng * sizeof( int );
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &aggr_cnt_array, nbytes, "AGI");
      for ( i = 0; i < aggr_cnt_leng; i++ ) aggr_cnt_array[i] = 0;
   } else
      aggr_cnt_array = NULL;

   /* ============================================================= */
   /* Phase 1  :                                                    */
   /*    for all nodes, form a new aggregate with its neighbors     */
   /*    if the number of its neighbors having been aggregated does */
   /*    not exceed a given threshold                               */
   /*    (max_neigh_selected = 0 ===> Vanek's scheme)               */
   /* ============================================================= */

   if ( ordering == 1 )
   {
      nbytes = (Nrows / num_PDE_eqns + 1) * sizeof(int);
      ML_memory_alloc((void**) &randomVector, nbytes, "AGJ");
      for (i = 0; i < Nrows/num_PDE_eqns; i++) randomVector[i] = i;
      ML_randomize(Nrows/num_PDE_eqns, randomVector);
   } 
   else if ( ordering == 2 )
   {
      new_node = (ML_Node *) malloc(sizeof(ML_Node));      
      new_node->node_id = 0;
      node_head = new_node;
      node_tail = new_node;
      new_node->next = NULL;
   }
   
   inode2 = 0;
   while ( inode2 < Nrows/num_PDE_eqns )
   {

      /* pick the next node to aggregate */

      if      ( ordering == 0 ) inode = inode2++;
      else if ( ordering == 1 ) inode = randomVector[inode2++];
      else if ( ordering == 2 ) 
      {
         if ( node_head == NULL ) 
         {
            for ( jnode = 0; jnode < Nrows/num_PDE_eqns; jnode++ ) 
            {
               if ( aggr_stat[jnode] == ML_AGGR_READY )
               { 
                  new_node = (ML_Node *) malloc(sizeof(ML_Node));      
                  new_node->node_id = jnode;
                  node_head = new_node;
                  node_tail = new_node;
                  new_node->next = NULL;
                  break;
               }
            }
         }
         if ( node_head == NULL ) break;
         new_node = node_head;
         inode = new_node->node_id;
         node_head = new_node->next;
         free(new_node);
      }

      /* consider further only if the node is in READY mode */

      if ( aggr_stat[inode] == ML_AGGR_READY ) 
      {
         length = amal_mat_indx[inode+1] - amal_mat_indx[inode] + 1;
         supernode = (ML_SuperNode *) malloc(sizeof(ML_SuperNode));      
         supernode->list = (int*) malloc(length*sizeof(int));

         if ((supernode->list) == NULL) 
         {
            printf("Error:couldn't allocate memory for supernode! %d\n",
                            length);
            exit(1);
         }

         supernode->maxlength = length;
         supernode->length = 1;
         supernode->list[0] = inode;
         select_flag = 1;

         /* count the number of neighbors that have been aggregated */

         count = 0;
         for (jnode=amal_mat_indx[inode];jnode<amal_mat_indx[inode+1];jnode++) 
         {
            index = amal_mat_indx[jnode];
            if ( index < Nrows/num_PDE_eqns ) 
            {
               if ( aggr_stat[index] == ML_AGGR_READY || 
                    aggr_stat[index] == ML_AGGR_NOTSEL ) 
               {
                  supernode->list[supernode->length++] = index;
               } 
               else if ( aggr_stat[index] != ML_AGGR_BDRY ) 
               {
                  count++;
               }
            }
         }

         /* if there are too many neighbors aggregated, don't do this one */

         if ( count > max_neigh_selected ) select_flag = 0;

         if (select_flag != 1 || supernode->length < min_nodes_per_aggregate) 
         {
            aggr_stat[inode] = ML_AGGR_NOTSEL;
            free( supernode->list );
            free( supernode );
            if ( ordering == 2 )
            {
               for (jnode=amal_mat_indx[inode];jnode<amal_mat_indx[inode+1];
                    jnode++)
               {
                  index = amal_mat_indx[jnode];
                  if ( aggr_stat[index] == ML_AGGR_READY )
                  {
                     new_node = (ML_Node *) malloc(sizeof(ML_Node));
                     new_node->node_id = index;
                     new_node->next = NULL;
                     if ( node_head == NULL )
                     {
                        node_head = new_node;
                        node_tail = new_node;
                     } else {
                        node_tail->next = new_node;
                        node_tail = new_node;
                     }
                  }
               }
            }
         } 
         else 
         {
            for ( j = 0; j < supernode->length; j++ ) 
            {
               jnode = supernode->list[j];
               aggr_stat[jnode] = ML_AGGR_SELECTED;
               aggr_index[jnode] = aggr_count;
               if ( ordering == 2 )
               {
                  for (kk=amal_mat_indx[jnode];kk<amal_mat_indx[jnode+1];kk++) 
                  {
                     if ( aggr_stat[amal_mat_indx[kk]] == ML_AGGR_READY )
                     { 
                        new_node = (ML_Node *) malloc(sizeof(ML_Node));      
                        new_node->node_id = amal_mat_indx[kk];
                        new_node->next = NULL;
                        if ( node_head == NULL )
                        {
                           node_head = new_node;
                           node_tail = new_node;
                        } else {
                           node_tail->next = new_node;
                           node_tail = new_node;
                        }
                     }
                  } 
               } 
            }
            supernode->next = NULL;
            supernode->index = aggr_count;
            if ( aggr_count == 0 ) 
            {
               aggr_head = supernode;
               aggr_curr = supernode;
            } 
            else 
            {
               aggr_curr->next = supernode;
               aggr_curr = supernode;
            } 
            aggr_cnt_array[aggr_count++] = supernode->length;
            if ( aggr_count >= aggr_cnt_leng ) 
            {
               itmp_array = aggr_cnt_array;
               aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
               nbytes = aggr_cnt_leng * sizeof( int );
               ML_memory_alloc((void**) &aggr_cnt_array, nbytes, "AGK");
               for ( k = 0; k < aggr_count; k++ )
                  aggr_cnt_array[k] = itmp_array[k];
               ML_memory_free((void**) &itmp_array);
            }
         }
      }
   }
   if ( ordering == 1 ) ML_memory_free((void**) &randomVector);
   else if ( ordering == 2 ) 
   {
      while ( node_head != NULL )
      {
         new_node = node_head;
         node_head = new_node->next;
         free( new_node );
      }
   }

   m = 0;
   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) 
      if ( aggr_stat[i] == ML_AGGR_READY ) m++;
   k = ML_Comm_GsumInt( comm, m);
   if ( k > 0 && mypid == 0 )
      printf("Aggregation(U) : Phase 1 (WARNING) - %d READY nodes left\n",k);
   m = 0;
   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) 
      if ( aggr_stat[i] == ML_AGGR_SELECTED ) m++;
   k = ML_Comm_GsumInt( comm, m);
   m = ML_Comm_GsumInt( comm, Nrows/num_PDE_eqns);
   j = ML_Comm_GsumInt( comm, aggr_count );
#ifdef ML_AGGR_PARTEST
	 global_ncoarse=m;
#endif
#ifdef ML_AGGR_DEBUG
   if ( mypid == 0 ) {
      printf("Aggregation(U) : Phase 1 - nodes aggregated = %d (%d)\n",k,m);
      printf("Aggregation(U) : Phase 1 - total aggregates = %d \n",j);
   }
#endif 
   /* ============================================================= */
   /* Phase 2 : aggregate the rest of the nodes into one of the     */
   /*           existing LOCAL aggregates. (attach_scheme)          */
   /* ============================================================= */
if (phase2_off == 0) {


   count = 0;
   for ( inode = 0; inode < Nrows/num_PDE_eqns; inode++ ) 
   {
      /* ---------------------------------------------------------- */
      /* for all nodes that have not been aggregated                */
      /* ---------------------------------------------------------- */

      if ( aggr_stat[inode] == ML_AGGR_NOTSEL || 
           aggr_stat[inode] == ML_AGGR_READY ) 
      {
         if ( attach_scheme == ML_AGGR_MINRANK ) 
         {
            /* search for a neighboring aggregate that has the fewest */
            /* number of nodes                                        */

            search_flag = 0;
            mincount = 100000;
            for (jnode=amal_mat_indx[inode]; jnode<amal_mat_indx[inode+1]; 
                 jnode++) 
            {
               index = amal_mat_indx[jnode];
               if ( index < Nrows/num_PDE_eqns ) 
               {
                  if ( aggr_stat[index] == ML_AGGR_SELECTED ) 
                  {
                     search_flag = 1;
                     m = aggr_index[index];
                     if ( aggr_cnt_array[m] < mincount ) 
                     {
                        mincount = aggr_cnt_array[m];
                        k = index;
                     }
                  }
               }
            }
            if ( search_flag == 1 ) 
            {
               index = k;
               m = aggr_index[index];
            }

         } 
         else if ( attach_scheme == ML_AGGR_MAXLINK ) 
         {
            /* search for a neighboring aggregate that has the most   */
            /* connection to my node                                  */

            search_flag = 0;
            length = amal_mat_indx[inode+1] - amal_mat_indx[inode];
            nbytes = length * sizeof( int );
            if ( nbytes > 0 )
               ML_memory_alloc((void**) &int_buf, nbytes, "AGR");
            length = 0; 
            for (jnode=amal_mat_indx[inode]; jnode<amal_mat_indx[inode+1]; 
                 jnode++) 
            {
               index = amal_mat_indx[jnode];
               if ( aggr_index[index] >= 0 ) 
                  int_buf[length++] = aggr_index[index];
            }
            ML_sort(length, int_buf);
            m = -1;
            maxcount = 0;
            if ( length > 0 ) {k = int_buf[0]; j = 1; m = k;}
            for ( jnode = 1; jnode < length; jnode++ ) 
            {
               if ( int_buf[jnode] == k ) j++;
               else 
               {
                  if ( j > maxcount ) 
                  {
                     maxcount = j;
                     m = k;
                  }
                  k = int_buf[jnode];
                  j = 1;
               }
            }
            if ( m >= 0 ) search_flag = 1;
            if ( nbytes > 0 ) ML_memory_free((void**) &int_buf);
         } else {
            printf("ML_Aggregate_CoarsenUncoupled error : invalid scheme.\n");
            exit(1);
         }

         /* if found, add the node to the existing aggregate       */

         if ( search_flag == 1 ) 
         { 
            aggr_cnt_array[m]++;
            aggr_index[inode] = m;
            aggr_stat[inode] = ML_AGGR_SELECTED2;
            /*###### don't need for now only
            supernode = aggr_head;
            for ( j = 0; j < m; j++ ) supernode = supernode->next;
            if ( supernode->length >= supernode->maxlength ) {
               length = supernode->length + 1;
               supernode->maxlength = length;
               itmp_array = supernode->list;
               supernode->list = (int*) malloc(length*sizeof(int));
               for ( j = 0; j < length-1; j++ ) 
                  supernode->list[j] = itmp_array[j];
               free( itmp_array );
            }
            supernode->list[supernode->length++] = inode;
            */
            count++;
         } 
         /*
         else {
            printf("Aggregation (U) : something wrong, node = %d\n", inode);
            printf("Aggregation (U) : please report to ML developers.\n");
         }
         */
      }
   }
}
   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) 
      if (aggr_stat[i] == ML_AGGR_SELECTED2) aggr_stat[i] = ML_AGGR_SELECTED;

   m = 0;
   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) 
      if ( aggr_stat[i] == ML_AGGR_SELECTED ) m++;
   k = ML_Comm_GsumInt( comm, m);
   m = ML_Comm_GsumInt( comm, Nrows/num_PDE_eqns);
   j = ML_Comm_GsumInt( comm, aggr_count );
#ifdef ML_AGGR_DEBUG
   if ( mypid == 0 ) {
      printf("Aggregation(U) : Phase 2 - nodes aggregated = %d (%d)\n",k,m);
      printf("Aggregation(U) : Phase 2 - total aggregates = %d \n",j);
   }
#endif

   /* ============================================================= */
   /* Phase 3 : for the un-aggregated nodes, form a new aggregate   */
   /* ============================================================= */

   for ( inode = 0; inode < Nrows/num_PDE_eqns; inode++ ) 
   {
      /* if not considered, and is a border node */

      if (aggr_stat[inode] == ML_AGGR_READY || 
          aggr_stat[inode] == ML_AGGR_NOTSEL ) 
      {
         /* examine all of its neighbors */

         count = 1;
         for (jnode = amal_mat_indx[inode]; jnode < amal_mat_indx[inode+1]; 
              jnode++) 
         {
            index = amal_mat_indx[jnode];
            if ( index < Nrows/num_PDE_eqns && 
                 aggr_stat[index] != ML_AGGR_SELECTED ) 
               count++;
         }
         length = amal_mat_indx[inode+1] - amal_mat_indx[inode];

         /* if enough neighbors have not been aggregated, form one */

         supernode = (ML_SuperNode *) malloc(sizeof(ML_SuperNode));      
         supernode->list = (int*) malloc(count*sizeof(int));
         if ((supernode->list) == NULL) 
         {
            printf("ML_Aggregate_Coarsen - couldn't allocate memory.\n");
            exit(1);
         }

         supernode->maxlength = count;
         supernode->length = 1;
         supernode->list[0] = inode;

         for (jnode = amal_mat_indx[inode]; jnode < amal_mat_indx[inode+1]; 
              jnode++) 
         {
            index = amal_mat_indx[jnode];
            if ( index < Nrows/num_PDE_eqns && 
                 aggr_stat[index] != ML_AGGR_SELECTED && 
                 aggr_stat[index] != ML_AGGR_BDRY ) 
               supernode->list[supernode->length++] = index;
         }
         for ( j = 0; j < supernode->length; j++ ) 
         {
            jnode = supernode->list[j];
            aggr_stat[jnode] = ML_AGGR_SELECTED;
            aggr_index[jnode] = aggr_count;
         }
         supernode->next = NULL;
         supernode->index = aggr_count;
         if ( aggr_count == 0 ) 
         {
            aggr_head = supernode;
            aggr_curr = supernode;
         } 
         else 
         {
            aggr_curr->next = supernode;
            aggr_curr = supernode;
         } 
         aggr_cnt_array[aggr_count++] = supernode->length;
         if ( aggr_count >= aggr_cnt_leng ) 
         {
            itmp_array = aggr_cnt_array;
            aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
            nbytes = aggr_cnt_leng * sizeof( int );
            ML_memory_alloc((void**) &aggr_cnt_array, nbytes, "AGL");
            for ( k = 0; k < aggr_count; k++ )
               aggr_cnt_array[k] = itmp_array[k];
            ML_memory_free((void**) &itmp_array);
         }
      }
   }

   m = 0;
   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) 
      if ( aggr_stat[i] == ML_AGGR_SELECTED ) m++;
   k = ML_Comm_GsumInt( comm, m);
   m = ML_Comm_GsumInt( comm, Nrows/num_PDE_eqns);
   j = ML_Comm_GsumInt( comm, aggr_count );

#ifdef ML_AGGR_DEBUG
   if ( mypid == 0 ) 
   {
      printf("Aggregation(U) : Phase 3 - nodes aggregated = %d (%d)\n",k,m);
      printf("Aggregation(U) : Phase 3 - total aggregates = %d \n",j);
   }
#endif

   /* ============================================================= */
   /* check for error                                               */
   /* ============================================================= */

   m = 0;
   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) 
      if (aggr_stat[i] != ML_AGGR_SELECTED && aggr_stat[i] != ML_AGGR_BDRY) m++;
   k = ML_Comm_GsumInt( comm, m);
   if ( k > 0 && mypid == 0 ) 
   {
      printf("Aggregation (U) error : not all nodes processed.\n");
      exit(1);
   }

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */

   Ncoarse = aggr_count * nullspace_dim;
 
   /* ------------------------------------------------------------- */
   /* check and copy aggr_index                                     */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) 
   {
      if ( (aggr_index[i] < 0 && aggr_stat[i] != ML_AGGR_BDRY) || 
            aggr_index[i] >= aggr_count) 
      {
         printf("Aggregation (U) error : index out of bound %d(%d)\n", 
                aggr_index[i], aggr_count);
      }
   }
   level = ml_ag->cur_level;
   nbytes = Nrows * sizeof( int );
   ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), nbytes, "AGv");
   count = aggr_count;
   for ( i = 0; i < Nrows/num_PDE_eqns; i++ ) 
   {
      if ( aggr_index[i] >= 0 )
      {
         for ( j = 0; j < num_PDE_eqns; j++ ) 
            ml_ag->aggr_info[level][i*num_PDE_eqns+j] = aggr_index[i];
      }
      else
      {
         for ( j = 0; j < num_PDE_eqns; j++ ) 
            ml_ag->aggr_info[level][i*num_PDE_eqns+j] = count;
         count++;
      }
   }
   /*ml_ag->aggr_count[level] = aggr_count; */
   ml_ag->aggr_count[level] = count; /* for relaxing boundary points */

#ifdef ML_AGGR_TEST
   if ( ml_ag->begin_level == ml_ag->cur_level && nprocs > 1 )
   {
      itmp_array  = (int *) malloc( nprocs * sizeof(int));
      itmp_array2 = (int *) malloc( nprocs * sizeof(int));
      for ( i = 0; i < nprocs; i++ ) itmp_array[i] = 0;
      itmp_array[mypid] = Nrows;
      ML_gsum_vec_int(itmp_array, itmp_array2, nprocs, comm);
      Nrows_offset = 0;
      for ( i = 0; i < mypid; i++ ) Nrows_offset += itmp_array[i];
      for ( i = 0; i < nprocs; i++ ) itmp_array[i] = 0;
      itmp_array[mypid] = count;
      ML_gsum_vec_int(itmp_array, itmp_array2, nprocs, comm);
      naggr_offset = 0;
      for ( i = 0; i < mypid; i++ ) naggr_offset += itmp_array[i];
      free(itmp_array);
      free(itmp_array2);
      sprintf(zfn, "mlaggr.out.%d",mypid);
      zfp = fopen(zfn,"w");
      fprintf(zfp, "aggr_count = %d \n", count);
      for ( i = 0; i < Nrows; i++ ) 
         fprintf(zfp, "%d \n", ml_ag->aggr_info[level][i]+naggr_offset);
      fclose(zfp);
      ML_gsum_vec_int(&i, &j, 1, comm);
      ML_gsum_vec_int(&i, &j, 1, comm);
      exit(1); 
   }
   if ( ml_ag->begin_level == ml_ag->cur_level && nprocs == 1 )
   {
      sprintf(zfn, "mlaggr.out");
      zfp = fopen(zfn,"r");
      for ( i = 0; i < Nrows; i++ ) 
         fscanf(zfp, "%d", &(ml_ag->aggr_info[level][i]));
      fclose(zfp);
      aggr_count = 0;
      for ( i = 0; i < Nrows; i++ ) 
         if ( ml_ag->aggr_info[level][i] > aggr_count )
            aggr_count = ml_ag->aggr_info[level][i];
      aggr_count++;
      Ncoarse = aggr_count * nullspace_dim;
      for ( i = 0; i < Nrows; i++ ) 
         aggr_index[i] = ml_ag->aggr_info[level][i];
      ml_ag->aggr_count[level] = aggr_count;
      printf("ML_Aggregate : total no. of aggregates input = %d\n",aggr_count);
   }
#endif

   /* ------------------------------------------------------------- */
   /* set up the space for storing the new operator and null space  */
   /* ------------------------------------------------------------- */

   nbytes = ( Nrows + 1 ) * sizeof(int); 
   ML_memory_alloc((void**)&(new_ia), nbytes, "AGM");
   nbytes = Nrows * nullspace_dim * sizeof(int);  
   ML_memory_alloc((void**)&(new_ja), nbytes, "AGN");
   nbytes = Nrows * nullspace_dim * sizeof(double); 
   ML_memory_alloc((void**)&(new_val), nbytes, "AGO");
   nbytes = Ncoarse * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&(new_null),nbytes,"AGX");
   for (i = 0; i < Ncoarse*nullspace_dim; i++) new_null[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* initialize the row pointer for the CSR prolongation operator  */
   /* (each roll will have at most nullspace_dim nonzero entries)   */
   /* ------------------------------------------------------------- */

   for (i = 0; i <  Nrows*nullspace_dim; i++) new_ja[i] = 0;
   for (i = 0; i <  Nrows*nullspace_dim; i++) new_val[i]= 0.0;
   for (i = 0; i <= Nrows; i++) new_ia[i] = i * nullspace_dim;
/* trying this when a Dirichlet row is taken out */
/*
j = 0;
for (i = 0; i <= Nrows; i++) {
   new_ia[i] = j;
   if (aggr_index[i] != -1) j += nullspace_dim;
}
*/

   /* ------------------------------------------------------------- */
   /* temporary variables for use in subsequent processing          */
   /* ------------------------------------------------------------- */

   nbytes = aggr_count * sizeof(int);
   ML_memory_alloc((void**)&agg_sizes,     nbytes,"AGi");
   ML_memory_alloc((void**)&agg_sizes_cum, nbytes, "AGj");

   /* ------------------------------------------------------------- */
   /* fill the temporary variables and also find the maximum        */
   /* aggregate size for allocation of qr_tmp                       */
   /* ------------------------------------------------------------- */

   for (i = 0; i < Nrows/num_PDE_eqns; i++) 
   {
      if (aggr_index[i] >= 0 && aggr_index[i] < aggr_count) 
         agg_sizes[aggr_index[i]] += num_PDE_eqns;
      else if (aggr_index[i] != -1) 
      {
         printf("%d : CoarsenUncoupled - wrong index %d(%d)\n",mypid,
                      aggr_index[i], aggr_count);
         exit(1);
      }
   }
   max_agg_size = agg_sizes[0];
   if ( aggr_count > 0 ) agg_sizes_cum[0] = 0;
   for (i = 0; i < aggr_count-1; i++)
   {
      agg_sizes_cum[i+1] = agg_sizes_cum[i] + agg_sizes[i];
      if (agg_sizes[i+1] > max_agg_size) max_agg_size = agg_sizes[i+1];
   }
   ML_memory_free((void**)&agg_sizes_cum);

   /* ------------------------------------------------------------- */
   /* generate an array to store which aggregate has which rows.Then*/
   /* loop through the rows of A checking which aggregate each row  */
   /* is in, and adding it to the appropriate spot in rows_in_aggs  */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**)&rows_in_aggs,aggr_count*sizeof(int*),"MLt");
   for (i = 0; i < aggr_count; i++) 
      rows_in_aggs[i] = (int *) malloc(agg_sizes[i]*sizeof(int));
   if (rows_in_aggs[aggr_count-1] == NULL) 
   {
      printf("Error: couldn't allocate memory in CoarsenUncoupled\n");
      exit(1);
   }
   for (i = 0; i < aggr_count; i++) agg_sizes[i] = 0;
   for (i = 0; i < Nrows/num_PDE_eqns; i++) 
   {
      if ( aggr_index[i] >= 0 )
      {
         for (j = 0; j < num_PDE_eqns; j++)
         {
            index = agg_sizes[aggr_index[i]]++; 
            rows_in_aggs[aggr_index[i]][index] = i*num_PDE_eqns + j;
         }
      }
   }

   /* ------------------------------------------------------------- */
   /* allocate work arrays for QR factorization                     */
   /* work and lwork are needed for lapack's QR routine.  These     */
   /* settings seemed easiest since I don't quite understand        */
   /* what they do, but may want to do something better here later  */
   /* ------------------------------------------------------------- */

   nbytes = max_agg_size * num_PDE_eqns * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&qr_tmp, nbytes, "AGU");
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&tmp_vect, nbytes, "AGT");

   lwork  = nullspace_dim;
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&work, nbytes, "AGk");

   for (i = 0; i < aggr_count; i++) 
   {
      /* set up the matrix we want to decompose into Q and R: */

      if (nullspace_vect == NULL) 
      {
         for (j = 0; j < agg_sizes[i]; j++)
         {
            for (k = 0; k < nullspace_dim; k++)
               if ( (rows_in_aggs[i][j]) % num_PDE_eqns == k)
                  qr_tmp[k*agg_sizes[i] + j] = 1.0;
               else
                  qr_tmp[k*agg_sizes[i] + j] = 0.0;
         }
      }
      else 
      {
         for (k = 0; k < nullspace_dim; k++)
            for (j = 0; j < agg_sizes[i]; j++)
               qr_tmp[k*agg_sizes[i] + j] = 
                  nullspace_vect[ k*Nrows + rows_in_aggs[i][j] ];
      }

      /* now calculate QR using an LAPACK routine */

      MLFORTRAN(dgeqrf)(&(agg_sizes[i]), &nullspace_dim, qr_tmp, 
                        &(agg_sizes[i]), tmp_vect, work, &lwork, &info);
      if (info != 0)
         pr_error("Error in CoarsenUncoupled : dgeqrf returned a non-zero\n");

      if (work[0] > lwork) 
      {
         lwork=(int) work[0]; 
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGk");
      }
      else lwork=work[0];
		 
      /* the upper triangle of qr_tmp is now R, so copy that into the 
         new nullspace */

      for (j = 0; j < nullspace_dim; j++)
         for (k = j; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse] = qr_tmp[j+agg_sizes[i]*k];
		 
      /* to get this block of P, need to run qr_tmp through another LAPACK 
         function: */

      MLFORTRAN(dorgqr)(&(agg_sizes[i]), &nullspace_dim, &nullspace_dim, 
              qr_tmp, &(agg_sizes[i]), tmp_vect, work, &lwork, &info);
      if (info != 0)
         pr_error("Error in CoarsenUncoupled: dorgqr returned a non-zero\n");

      if (work[0] > lwork) 
      {
         lwork=(int) work[0]; 
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGk");
      }
      else lwork=work[0];
		 
      /* now copy Q over into the appropriate part of P: */
      /* The rows of P get calculated out of order, so I assume the Q is 
         totally dense and use what I know of how big each Q will be to 
         determine where in ia, ja, etc each nonzero in Q belongs.  If I 
         did not assume this, I would have to keep all of P in memory in 
         order to determine where each entry should go */

      for (k = 0; k < nullspace_dim; k++) 
      {
         largest = 0.; thesign = 1.;
         for (j = 0; j < agg_sizes[i]; j++)
         {
            if ( dabs(qr_tmp[ k*agg_sizes[i] + j ]) > largest ) {
               largest = dabs(qr_tmp[ k*agg_sizes[i] + j ]);
               if ( qr_tmp[ k*agg_sizes[i] + j ] < 0.0) thesign = -1.;
               else thesign = 1.;
            }
         }
         for (j = 0; j < agg_sizes[i]; j++)
         {
            index = new_ia[rows_in_aggs[i][j]] + k;
            new_ja [index] = i * nullspace_dim + k;
            new_val[index] = thesign*qr_tmp[ k*agg_sizes[i] + j ];
         }
      }
   }
	 
   /*for (j = 0; j < Ncoarse; j++)
    *{
    *   printf("null %4d : ", j);
    *   for (k = 0; k < nullspace_dim; k++)
    *      printf("(%4d %13.7e) ", k, new_null[k*Ncoarse+j]);
    *   printf("\n");
    *} */

   ML_Aggregate_Set_NullSpace(ml_ag, num_PDE_eqns, nullspace_dim, 
                              new_null, Ncoarse);
   ML_memory_free( (void **) &new_null);

   /* ------------------------------------------------------------- */
   /* compress the prolongation operator                            */
   /* ------------------------------------------------------------- */

   kk    = new_ia[0];
   index = kk;
   for (i = 0; i < Nrows; i++)
   {
      for (j = kk; j < new_ia[i+1]; j++ )
      {
         if ( new_val[j] != 0.0 )
         {
            new_val[index]  = new_val[j];  
            new_ja[index++] = new_ja[j];  
         }
      }
      if ( index == new_ia[i] ) 
      {
         new_val[index] = new_val[kk]; new_ja[index++] = new_ja[kk];
      }
      kk = new_ia[i+1];
      new_ia[i+1] = index;
   }
   /*if ( index > 2 * Nrows )
    *{
    *   int_buf  = new_ja;
    *   dble_buf = new_val;
    *   nbytes = index * sizeof(int);  
    *   ML_memory_alloc((void**)&(new_ja), nbytes, "AGN");
    *   nbytes = index * nullspace_dim * sizeof(double); 
    *   ML_memory_alloc((void**)&(new_val), nbytes, "AGO");
    *   for ( i = 0; i < index; i++ )
    *   {
    *      new_ja[i] = int_buf[i];
    *      new_val[i] = dble_buf[i];
    *   }
    *   ML_memory_free((void**)&(int_buf));
    *   ML_memory_free((void**)&(dble_buf));
    *} 
    *for (j = 0; j < Nrows; j++)
    *   for (i = new_ia[j]; i < new_ia[j+1]; i++)
    *      printf("P(%4d,%4d) = %e\n", j, new_ja[i], new_val[i]);
   */

#if defined(MATLAB) && defined(PARTEST)
   proc=comm->ML_mypid;
	 
   sprintf(pfn,"PmatrixProc%dNc%d.dat", proc, global_ncoarse);
   sprintf(nfn,"NullSpaceProc%dNc%d.dat", proc, global_ncoarse);
   sprintf(afn,"AggregatesProc%dNc%d.dat", proc, global_ncoarse);
	 
   pfp=fopen(pfn,"w");
   nfp=fopen(nfn,"w");
   afp=fopen(afn,"w");
	 
   for (i=0; i<Nrows; i++) {
      globrow=local_to_global_row(i);
      if (globrow<0)
         printf("Argh!, index=%d, proc=%d, globrow=%d\n", i, proc, globrow); 
      for (j=0; j < nullspace_dim; j++) {
         fprintf(pfp,"%d %d %f %d\n", globrow, new_ja[i*nullspace_dim+j], 
         new_val[i*nullspace_dim+j], i);
      }
      fprintf(afp, "%d %d %d\n", aggr_index[i/num_PDE_eqns], globrow, i);
   }

   for (i = 0; i < Ncoarse; i++)
   for (j = 0; j < nullspace_dim; j++)
      fprintf(nfp,"%d %d %f %d\n", i, j, new_null[j*Ncoarse+i], i);
	 
   fclose(pfp);
   fclose(nfp);
   fclose(afp);
#endif

#ifdef MATLAB2
   sprintf(afn,"AggregatesNc%d.dat", Ncoarse);
   afp=fopen(afn,"w");
   for (i=0; i<Nrows; i++) {
      fprintf(afp, "%d\n", aggr_index[i/num_PDE_eqns]);
   }
   fclose(afp);
#endif

   ML_memory_alloc((void**) &csr_data, sizeof(struct ML_CSR_MSRdata), "AGP");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;
   (*Pmatrix) = ML_Operator_Create(comm);
   ML_Operator_Set_ApplyFuncData( *Pmatrix, Ncoarse, Nrows, ML_EMPTY,
                                  csr_data, Nrows, NULL, 0);
   (*Pmatrix)->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   ML_memory_alloc((void**) &aggr_comm, sizeof(ML_Aggregate_Comm), "AGQ");
   aggr_comm->comm = comm;
   aggr_comm->N_send_neighbors = 0;
   aggr_comm->N_recv_neighbors = 0;
   aggr_comm->send_neighbors = NULL;
   aggr_comm->recv_neighbors = NULL;
   aggr_comm->send_leng = NULL;
   aggr_comm->recv_leng = NULL;
   aggr_comm->send_list = NULL;
   aggr_comm->local_nrows = Ncoarse;
   
   ML_CommInfoOP_Generate( &((*Pmatrix)->getrow->pre_comm), 
                           ML_Aggregate_ExchangeBdry, aggr_comm, 
                           comm, Ncoarse, m);
   ML_Operator_Set_Getrow((*Pmatrix), ML_EXTERNAL, Nrows, CSR_getrows);
   ML_Operator_Set_ApplyFunc((*Pmatrix), ML_INTERNAL, CSR_matvec);

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_memory_free((void**) &mat_indx);
   if (num_PDE_eqns !=1) ML_memory_free((void**) &amal_mat_indx);
   ML_memory_free((void**) &aggr_index);
   ML_memory_free((void**) &aggr_stat);
   ML_memory_free((void**) &aggr_cnt_array);
   ML_memory_free((void**) &aggr_comm);
   aggr_curr = aggr_head;
   while ( aggr_curr != NULL ) 
   {
      supernode = aggr_curr;
      aggr_curr = aggr_curr->next;
      if ( supernode->maxlength > 0 ) free( supernode->list );
      free( supernode );
   }
   ML_memory_free((void**)&agg_sizes);
   for (i = 0; i < aggr_count; i++) free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   ML_memory_free((void**)&qr_tmp);
   ML_memory_free((void**)&tmp_vect);
   ML_memory_free((void**)&work);

   return Ncoarse;
}

/* ******************************************************************** */
/* ******************************************************************** */
/* The following subroutine is useful to ensure coupling between        */
/* processors.  However, it has not been updated as good as the         */
/* ML_Aggregate_CoarsenUncoupled subroutine.                            */
/*               --- Charles Tong, December 29, 1999                    */
/* ******************************************************************** */
/* ******************************************************************** */
/* ******************************************************************** */
/* construct the tentative prolongator allowing aggregate to cross      */
/* processor boundaries                                                 */
/* -------------------------------------------------------------------- */

int ML_Aggregate_CoarsenCoupled( ML_Aggregate *ml_ag, ML_Operator *Amatrix, 
                                 ML_Operator **Pmatrix, ML_Comm *comm)
{
   int     i, j = 0, k, kk, m, inode, inode2, jnode, nbytes, length = 0;
   int     *mat_indx = NULL, jj, Nrows, exp_Nrows, *int_deg_list = NULL;
   int     *ext_deg_list = NULL, N_neighbors, *neighbors = NULL;
   int     *recv_leng = NULL, *send_leng = NULL, *send_list = NULL;
   int     total_recv_leng, total_send_leng, offset, msgtype;
   int     select_flag, aggr_count, index, mypid, loop_flag, mdiff;
   int     *aggr_index = NULL, search_flag, *itmp_array = NULL;
   int     *sendlist_proc = NULL, Ncoarse, count, mincount, *int_buf=NULL;
   int     *int_buf2, old_aggr_count, *aggr_stat = NULL, procnum, index2;
   int     nz_cnt, index4, oldleng, *new_send_leng = NULL, new_N_send;
   int     *new_send_neighbors = NULL, status, *new_send_list = NULL;
   int     max_count, *new_ia = NULL, *new_ja = NULL, new_Nrows;
   int     *new_recv_leng = NULL, exp_Ncoarse, new_N_recv, maxcount;
   int     *new_recv_neighbors = NULL, aggr_cnt_leng;
   int     *aggr_cnt_array = NULL, *randomVector = NULL, ordering;
   int     maxnnz_per_row=500, diff_level;
   int     attach_scheme, num_PDE_eqns, nullspace_dim;
   int     level, index3, count3, *recv_list = NULL, max_agg_size;
   int     **rows_in_aggs = NULL, lwork, info, *trackbc, bdry_flag, nn2;
   double  *diagonal = NULL, dcompare1;
   double  dcompare2, *new_val = NULL, epsilon;
   double  *dble_buf = NULL, *nullspace_vect = NULL, *qr_tmp = NULL;
   double  *tmp_vect = NULL, *work = NULL, *new_null=NULL, *comm_val=NULL;
   double  *dble_buf2;
   int     (*getrowfunc)(void *,int,int*,int,int*,double*,int*);
   ML_SuperNode          *aggr_head, *aggr_curr = NULL, *supernode;
   ML_Node               *node_head = NULL, *new_node;
   struct ML_CSR_MSRdata *csr_data;
   ML_Aggregate_Comm     *aggr_comm;
   ML_GetrowFunc         *getrow_obj;
   ML_CommInfoOP         *getrow_comm;
   USR_REQ               *request = NULL;

   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   mypid          = comm->ML_mypid;
   epsilon        = ml_ag->threshold;
   ordering       = ml_ag->ordering;
   attach_scheme  = ml_ag->attach_scheme;
   num_PDE_eqns   = ml_ag->num_PDE_eqns;
   nullspace_dim  = ml_ag->nullspace_dim;
   nullspace_vect = ml_ag->nullspace_vect;
   Nrows          = Amatrix->outvec_leng;

   /* ============================================================= */
   /* check the system size versus null dimension size              */
   /* ============================================================= */

   if ( Nrows % num_PDE_eqns != 0 )
   {
      printf("ML_Aggregate_CoarsenCoupled ERROR : Nrows must be multiples");
      printf(" of num_PDE_eqns.\n");
      exit(1);
   }
   diff_level = ml_ag->max_levels - ml_ag->cur_level - 1;
   if ( diff_level > 0 ) num_PDE_eqns = nullspace_dim; /* ## 12/20/99 */

   /* ============================================================= */
   /* set up the threshold for weight-based coarsening              */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level < 0 ) diff_level = - diff_level;
   if ( diff_level > 0 ) num_PDE_eqns = nullspace_dim; /* ## 6/14/00 */
   while ( diff_level > 0 ) {epsilon *= 0.5; diff_level--;}

#ifdef ML_AGGR_DEBUG
   if ( mypid == 0 )
   {
      printf("ML_Aggregate_CoarsenCoupled : current level = %d\n",
                                            ml_ag->cur_level);
      printf("ML_Aggregate_CoarsenCoupled : current eps = %e\n",epsilon);
   }
#endif
   epsilon = epsilon * epsilon;

   /* ============================================================= */
   /* fetch the getrow function for the incoming matrix             */
   /* ============================================================= */

   getrow_obj = Amatrix->getrow;
   if ( getrow_obj->ML_id == ML_EXTERNAL) 
   {
      getrowfunc = getrow_obj->external;
   } 
   else 
   {
      getrowfunc = getrow_obj->internal;
   }
   if ( getrowfunc == NULL ) 
   {
      printf("ML_Aggregate_CoarsenCoupled Error : getrowfunc is null.\n");
      exit(-1);
   }

   /* ============================================================= */
   /* allocate initial temporary storage space for getrow           */
   /* also allocate space for storing the diagonal                  */
   /* ============================================================= */

   nbytes = maxnnz_per_row * sizeof( int );
   ML_memory_alloc((void**) &int_buf, nbytes, "AGA");
   nbytes = maxnnz_per_row * sizeof( double );
   ML_memory_alloc((void**) &dble_buf, nbytes, "AGB");
   if ( Nrows > 0 )
   {
      nbytes = Nrows * sizeof( double );
      ML_memory_alloc((void**) &diagonal, nbytes, "AGC");
   }
   else diagonal = NULL;

   /* ============================================================= */
   /* fill in the diagonal array, also find out about the size of   */
   /* the incoming matrix (for allocation purpose)                  */
   /* ============================================================= */

   exp_Nrows = Nrows - 1;
   count = 0;
   for ( i = 0; i < Nrows; i++ ) 
   {
      diagonal[i]     = 0.0;
      while (getrowfunc(Amatrix->data,1,&i,maxnnz_per_row,int_buf, 
                        dble_buf, &m) == 0 ) 
      {
         ML_memory_free((void**) &int_buf);
         ML_memory_free((void**) &dble_buf);
         maxnnz_per_row = maxnnz_per_row * 2 + 1; 
         nbytes = maxnnz_per_row * sizeof( int );
         ML_memory_alloc((void**) &int_buf, nbytes, "AGD");
         nbytes = maxnnz_per_row * sizeof( double );
         ML_memory_alloc((void**) &dble_buf,  nbytes, "AGE");
      }
      for ( j = 0; j < m; j++ ) 
      {
         if ( int_buf[j] > exp_Nrows ) exp_Nrows = int_buf[j];
         if ( int_buf[j] == i )        diagonal[i] = dble_buf[j];
      }
      count += m;
      if ( diagonal[i] == 0.0 ) 
      {
         printf("%d : CoarsenCoupled WARNING - diag %d is 0.\n",mypid,i);
         count++;
      }
   }
   if ( int_buf  != NULL ) ML_memory_free((void**) &int_buf);
   if ( dble_buf != NULL ) ML_memory_free((void**) &dble_buf);
   exp_Nrows++;
   if ( epsilon == 0.0 && diagonal != NULL )
   {
      ML_memory_free((void**) &diagonal);
      diagonal = NULL;
   } 
   else 
   {
      dble_buf = diagonal;
      nbytes = exp_Nrows * sizeof(double);
      ML_memory_alloc((void**) &diagonal, nbytes, "AGF");
      for ( i = 0; i < Nrows; i++ ) diagonal[i] = dble_buf[i];
      for ( i = Nrows; i < exp_Nrows; i++ ) diagonal[i] = 0.0;
      ML_memory_free((void**) &dble_buf);
   }
#ifdef ML_AGGR_OBSOLETE
      /* The goal of this code is to set the diagonal such that */
      /* only the first 30 nonzeros are considered when doing   */
      /* aggregation.  I'm not 100% sure that it works ... but  */
      /* would like to keep this.                               */
      nbytes = maxnnz_per_row * sizeof( double );
      ML_memory_alloc((void**) &dble_buf,  nbytes, "AGy");
      for ( i = 0; i < Nrows; i++ ) {
         diagonal[i] = 1.0e-30;
         if ( mat_indx[i+1] - mat_indx[i] > 30) {
            for (j = mat_indx[i] ; j < mat_indx[i+1]; j++) {
               dble_buf[j-mat_indx[i]] = dabs(mat_val[j]);
            }
            ML_az_dsort( dble_buf, mat_indx[i+1]-mat_indx[i]);
            diagonal[i] =dble_buf[mat_indx[i+1]-mat_indx[i] - 30]/sqrt(epsilon);
         }
      }
      ML_memory_free((void**) &dble_buf);
#endif

   /* ============================================================= */
   /* build the diagonals of the expanded rows if epsilon != 0      */
   /* (diagonal elements are needed for pruning weak edges)         */
   /* ============================================================= */

   if ( epsilon != 0.0 )
   {
      getrow_comm = getrow_obj->pre_comm;
      if ( getrow_comm != NULL )
         ML_exchange_bdry(diagonal,getrow_comm,Nrows,comm,ML_OVERWRITE);
   }
      
   /* ============================================================= */
   /* construct 2 arrays : one to store the number of external      */
   /* links, and one to store the number of internal links for each */
   /* node.                                                         */
   /* ============================================================= */

   nbytes = Nrows * sizeof( int );
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &int_deg_list, nbytes, "AGG");
      ML_memory_alloc((void**) &ext_deg_list, nbytes, "AGH");
      ML_memory_alloc((void**) &trackbc,      nbytes, "AGI");
   } 
   else
   {
      int_deg_list = NULL;
      ext_deg_list = NULL;
      trackbc      = NULL;
   }

   /* ============================================================= */
   /* allocate temporary storage space for getrow                   */
   /* ============================================================= */

   nbytes = (count + 1) * sizeof( int );
   ML_memory_alloc((void**) &mat_indx, nbytes, "AGI");
   k = ML_Comm_GsumInt( comm, Nrows);
   m = ML_Comm_GsumInt( comm, count);
   nbytes = maxnnz_per_row * sizeof(int);
   ML_memory_alloc((void**) &int_buf, nbytes, "AGJ");
   nbytes = maxnnz_per_row * sizeof(double);
   ML_memory_alloc((void**) &dble_buf, nbytes, "AGK");

   k = ML_Comm_GsumInt( comm, Nrows);
   m = ML_Comm_GsumInt( comm, count);
#ifdef ML_AGGR_DEBUG
   if ( mypid == 0 )
      printf("Aggregation(C) : Total nnz = %d (Nrows=%d)\n",m,k);
#endif
   if ( ml_ag->operator_complexity == 0.0 )
   {
      ml_ag->fine_complexity = 1.0 * m;
      ml_ag->operator_complexity = 1.0 * m;
   }
   else
   {
      ml_ag->operator_complexity += 1.0 * m;
   }

   /* ============================================================= */
   /* extract the matrix using the getrow function                  */
   /* ============================================================= */

   count = Nrows + 1;
   mat_indx[0] = count;
   for ( i = 0; i < Nrows; i++ ) 
   {
      getrowfunc(Amatrix->data,1,&i,maxnnz_per_row,int_buf,dble_buf,&m);
      if ( m == 0 ) trackbc[i] = 1; else trackbc[i] = 0;
      nz_cnt = 0;
      int_deg_list[i] = 1;
      ext_deg_list[i] = 0;
      index3 = i / num_PDE_eqns * num_PDE_eqns;
      for (j = 0; j < m; j++) 
      {
         jnode = int_buf[j];
         if (((jnode - index3) < 0 || (jnode-index3) >= num_PDE_eqns) && 
               epsilon > 0.0)
         {
            dcompare1 = dble_buf[j] * dble_buf[j];
            if ( dcompare1 > 0.0 )
            {
               nz_cnt++;
               dcompare2 = dabs((diagonal[i] * diagonal[jnode]));
               if ( dcompare1 >= epsilon * dcompare2 )
                  mat_indx[count++] = int_buf[j];
            }
         } 
         else if (((jnode - index3) < 0 || (jnode-index3) >= num_PDE_eqns) && 
                  dble_buf[j] != 0.0)
         {
            mat_indx[count++] = int_buf[j];
            nz_cnt++;
         }
      }
      for (j = mat_indx[i]; j < count; j++) 
      {
         if ( mat_indx[j] >= Nrows ) ext_deg_list[i]++;
         else                        int_deg_list[i]++;
      }
      mat_indx[i+1] = count;
      ML_sort(mat_indx[i+1]-mat_indx[i], mat_indx+mat_indx[i]);
   }
   if ( int_buf  != NULL ) ML_memory_free((void**) &int_buf);
   if ( dble_buf != NULL ) ML_memory_free((void**) &dble_buf);
   if ( diagonal != NULL ) ML_memory_free((void**) &diagonal);

   /* ============================================================= */
   /* extract the communication information for the incoming matrix */
   /*  - N_neighbors, neighbors, send_leng, recv_leng, send_list    */
   /*  - total_send_leng, total_recv_leng                           */
   /* This is needed to establish aggregating order for border nodes*/
   /* ============================================================= */

   if ( num_PDE_eqns == 1 )
   {
      N_neighbors = getrow_obj->pre_comm->N_neighbors;
      nbytes = N_neighbors * sizeof( int );
      if ( nbytes > 0 ) 
      {
         ML_memory_alloc((void**) &neighbors,  nbytes, "AGL");
         ML_memory_alloc((void**) &recv_leng,  nbytes, "AGM");
         ML_memory_alloc((void**) &send_leng,  nbytes, "AGN");
      } 
      else 
      {
         neighbors = recv_leng = send_leng = NULL;
      }
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         neighbors[i] = getrow_obj->pre_comm->neighbors[i].ML_id;
         recv_leng[i] = getrow_obj->pre_comm->neighbors[i].N_rcv;
         send_leng[i] = getrow_obj->pre_comm->neighbors[i].N_send;
      }
      total_recv_leng = total_send_leng = 0;
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         total_recv_leng += recv_leng[i];
         total_send_leng += send_leng[i];
      }
      nbytes = total_send_leng * sizeof( int );
      if ( nbytes > 0 ) ML_memory_alloc((void**) &send_list,nbytes,"AGO");
      else              send_list = NULL;
      if ( total_recv_leng+Nrows != exp_Nrows ) 
      {
         printf("%d : ML_Aggregate_CoarsenCoupled - internal error.\n",mypid);
         printf("     lengths = %d %d \n",total_recv_leng+Nrows,exp_Nrows);
         exit(-1);
      }
      count = 0;
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         for (j = 0; j < send_leng[i]; j++)
            send_list[count++] = 
               getrow_obj->pre_comm->neighbors[i].send_list[j];
      }
      if ( count > total_send_leng )
      {
         printf("%d : CoarsenCoupled ERROR : count < total_send_leng\n",mypid);
         exit(1);
      }
   }
   else
   {
      /* ---------------------------------------------------------- */
      /* allocate storage for the communication information         */
      /* ---------------------------------------------------------- */

      N_neighbors = getrow_obj->pre_comm->N_neighbors;
      nbytes = N_neighbors * sizeof( int );
      if ( nbytes > 0 ) 
      {
         ML_memory_alloc((void**) &neighbors,  nbytes, "AGL");
         ML_memory_alloc((void**) &recv_leng,  nbytes, "AGM");
         ML_memory_alloc((void**) &send_leng,  nbytes, "AGN");
      } 
      else 
      {
         neighbors = recv_leng = send_leng = NULL;
      }
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         neighbors[i] = getrow_obj->pre_comm->neighbors[i].ML_id;
         recv_leng[i] = getrow_obj->pre_comm->neighbors[i].N_rcv;
         send_leng[i] = getrow_obj->pre_comm->neighbors[i].N_send;
      }
      total_recv_leng = 0;
      for ( i = 0; i < N_neighbors; i++ ) total_recv_leng += recv_leng[i];
      total_send_leng = 0;
      for ( i = 0; i < N_neighbors; i++ ) total_send_leng += send_leng[i];
      nbytes = total_send_leng * num_PDE_eqns * sizeof( int );
      if ( nbytes > 0 ) ML_memory_alloc((void**) &send_list,nbytes,"AGO");
      else              send_list = NULL;
      nbytes = total_recv_leng * sizeof( int );
      if ( nbytes > 0 ) ML_memory_alloc((void**) &recv_list,nbytes,"AGP");
      else              recv_list = NULL;

      /* ---------------------------------------------------------- */
      /* set up true external indices to be shipped to receive      */
      /* processors (true in view of that num_PDE_eqns can be > 1)  */
      /* ---------------------------------------------------------- */

      nbytes = Nrows * sizeof( int );
      if ( nbytes > 0 ) ML_memory_alloc((void**) &itmp_array,nbytes,"AGQ");
      count = 0;
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         for ( j = 0; j < Nrows; j++ ) itmp_array[j] = -1;
         count3 = 0;
         for (j = 0; j < send_leng[i]; j++)
         {
            index3 = getrow_obj->pre_comm->neighbors[i].send_list[j];
            index3 = index3 / num_PDE_eqns * num_PDE_eqns;
            for (k = 0; k < num_PDE_eqns; k++)
            {
               if ( itmp_array[index3+k] < 0 ) 
                  itmp_array[index3+k] = count3++;
            }
         }
         for (j = 0; j < send_leng[i]; j++)
         {
            send_list[count+j] = 
               getrow_obj->pre_comm->neighbors[i].send_list[j];
         }
         for ( j = 0; j < send_leng[i]; j++ ) 
         {
            index = send_list[count+j];
            if (itmp_array[index] >= 0) send_list[count+j] = itmp_array[index];
         }
         count += send_leng[i];
      }
      ML_memory_free((void**) &itmp_array);

      /* ---------------------------------------------------------- */
      /* send the adjusted indices to the receive processors        */
      /* ---------------------------------------------------------- */

      if ( N_neighbors > 0 )
         request = (USR_REQ *) malloc(N_neighbors*sizeof(USR_REQ));

      offset = 0;
      for (i = 0; i < N_neighbors; i++) 
      {
         msgtype = 2000;
         length = recv_leng[i] * sizeof( int );
         procnum = neighbors[i];
         comm->USR_irecvbytes((void *) &(recv_list[offset]),length,&procnum,
                              &msgtype, comm->USR_comm, request+i);
         offset += recv_leng[i];
      }
      offset = 0;
      for (i = 0; i < N_neighbors; i++) 
      {
         msgtype = 2000;
         length = send_leng[i] * sizeof( int );
         procnum = neighbors[i];
         comm->USR_sendbytes((void *) &(send_list[offset]),length,procnum, 
                              msgtype, comm->USR_comm);
         offset += send_leng[i];
      }
      offset = 0;
      for (i = 0; i < N_neighbors; i++) 
      {
         msgtype = 2000;
         length = recv_leng[i] * sizeof( int );
         procnum = neighbors[i];
         comm->USR_waitbytes((void *) &(recv_list[offset]),length,&procnum,
                              &msgtype, comm->USR_comm, request+i);
         for (j = 0; j < recv_leng[i]; j++) recv_list[offset+j] += offset; 
         offset += recv_leng[i];
      }
      if ( N_neighbors > 0 ) free( request );

      /* ---------------------------------------------------------- */
      /* adjust the matrix to reflect the index remapping           */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < Nrows; i++ )
      {
         for ( j = mat_indx[i]; j < mat_indx[i+1]; j++ )
            if ( mat_indx[j] >= Nrows ) 
               mat_indx[j] = recv_list[mat_indx[j]-Nrows] + Nrows;
      }
      ML_memory_free((void**) &recv_list);

      /* ---------------------------------------------------------- */
      /* update the send_list and send_leng's in line with remap    */
      /* ---------------------------------------------------------- */

      nbytes = Nrows * sizeof( int );
      if (nbytes > 0) ML_memory_alloc((void**) &itmp_array,nbytes,"AGR");
      total_send_leng = 0;
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         count = 0;
         for ( j = 0; j < Nrows; j++ ) itmp_array[j] = -1;
         for (j = 0; j < send_leng[i]; j++)
         {
            index3 = getrow_obj->pre_comm->neighbors[i].send_list[j];
            index3 = index3 / num_PDE_eqns * num_PDE_eqns;
            for (k = 0; k < num_PDE_eqns; k++)
               itmp_array[index3+k] = 0;
         }
         for ( j = 0; j < Nrows; j++ ) 
         {
            if ( itmp_array[j] == 0 ) send_list[total_send_leng+count++] = j;
         }
         send_leng[i] = count;
         total_send_leng += count;
      }
      total_send_leng = 0;
      for ( i = 0; i < N_neighbors; i++ ) total_send_leng += send_leng[i];
    
      ML_memory_free((void**) &itmp_array);

      /* ---------------------------------------------------------- */
      /* update other processors with the new communication pattern */
      /* ---------------------------------------------------------- */
     
      if ( N_neighbors > 0 )
         request = (USR_REQ *) malloc(N_neighbors*sizeof(USR_REQ));

      for (i = 0; i < N_neighbors; i++) 
      {
         msgtype = 2002;
         length = sizeof( int );
         procnum = neighbors[i];
         comm->USR_irecvbytes((void *) &(recv_leng[i]), length, &procnum,
                              &msgtype, comm->USR_comm, request+i);
      }
      for (i = 0; i < N_neighbors; i++) 
      {
         msgtype = 2002;
         length = sizeof( int );
         procnum = neighbors[i];
         comm->USR_sendbytes((void *) &(send_leng[i]), length, procnum, 
                              msgtype, comm->USR_comm);
      }
      for (i = 0; i < N_neighbors; i++) 
      {
         msgtype = 2002;
         length = sizeof( int );
         procnum = neighbors[i];
         comm->USR_waitbytes((void *) &(recv_leng[i]), length, &procnum,
                              &msgtype, comm->USR_comm, request+i);
      }
      if ( N_neighbors > 0 ) free( request );

      total_recv_leng = 0;
      for (i = 0; i < N_neighbors; i++) total_recv_leng += recv_leng[i];
      exp_Nrows = Nrows + total_recv_leng;;
   }

   /* ============================================================= */
   /* sendlist_proc is used to find out, in the aggregation process,*/
   /* which processor holds neighbors of my local nodes             */
   /* ============================================================= */

   nbytes = (N_neighbors + 1) * sizeof(int);
   ML_memory_alloc((void**) &sendlist_proc, nbytes, "AGS");
   sendlist_proc[0] = 0;
   for ( i = 1; i <= N_neighbors; i++ ) 
      sendlist_proc[i] = sendlist_proc[i-1] + recv_leng[i-1]; 

   /* ============================================================= */
   /* set up bookkeeping mechanism for aggregation                  */
   /* ============================================================= */

   aggr_count = 0;
   aggr_head = NULL;
   aggr_cnt_leng = Nrows / 5 + 2;
   nbytes = aggr_cnt_leng * sizeof( int );
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &aggr_cnt_array, nbytes, "AGT");
      for ( i = 0; i < aggr_cnt_leng; i++ ) aggr_cnt_array[i] = 0;
   } else
      aggr_cnt_array = NULL;

   /* ============================================================= */
   /* Construct an initial status array to store whether the nodes  */
   /* have been aggregated. aggr_index stores the aggregate number  */
   /* where this node has been aggregated into.                     */
   /* ============================================================= */

   nbytes = exp_Nrows * sizeof( int );
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &aggr_index, nbytes, "AGU");
      ML_memory_alloc((void**) &aggr_stat,  nbytes, "AGV");
   } else aggr_index = aggr_stat = NULL;
   for ( i = 0; i < Nrows; i++ ) 
   {
      if (trackbc[i] == 1) aggr_stat[i] = ML_AGGR_BDRY;
      else                 aggr_stat[i] = ML_AGGR_READY;
   }
   ML_memory_free((void **) &trackbc);
   for ( i = 0; i < Nrows; i+= num_PDE_eqns ) 
   {
      count = 0;
      for ( j = 0; j < num_PDE_eqns; j++ )
      {
         if ( aggr_stat[i] == ML_AGGR_BDRY ) count++;
      }
      if ( count != num_PDE_eqns )
         for (j = 0; j < num_PDE_eqns; j++) aggr_stat[j] = ML_AGGR_READY;
   } 
  
   for ( i = 0; i < exp_Nrows; i++ ) aggr_index[i] = -1;
   for ( i = Nrows; i < exp_Nrows; i++ ) aggr_stat[i] = 0;

   /* ============================================================= */
   /* allocate communication buffers for exchanging status info     */
   /* between processors during the aggregation step                */
   /* ============================================================= */

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf, nbytes, "AGW");
   else              int_buf = NULL;
   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf2, nbytes, "AGX");
   else              int_buf2 = NULL;

   /* ============================================================= */
   /* label all nodes that have to wait for its neighbors           */
   /* (If any of my neighbors reside in processors with processor   */
   /*  ID lower than my processor ID, then my node has to wait)     */
   /* ============================================================= */

   for ( inode = 0; inode < Nrows; inode++ ) 
   {
      if (ext_deg_list[inode] != 0) /* label for border nodes only */
      {
         for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1]; jnode++) 
         {
            index = mat_indx[jnode];
            mdiff = index - Nrows;

            /* ---------------------------------------------------- */
            /* search for the processor the node is coming from     */
            /* ---------------------------------------------------- */

            for ( k = 0; k <= N_neighbors; k++ )
               if ( mdiff < sendlist_proc[k] ) break;

            /* ---------------------------------------------------- */
            /* if the processor number < mypid, tag it with the     */
            /* neighbor processor with the smallest rank            */
            /* ---------------------------------------------------- */

            if ( k != 0 && neighbors[k-1] < mypid )
            {
               if ( aggr_stat[inode] < 0 )
                  aggr_stat[inode] = neighbors[k-1];
               else if ( neighbors[k-1] < aggr_stat[inode] )
                  aggr_stat[inode] = neighbors[k-1];
            }
         }
      }
   }

   /* ------------------------------------------------------------- */
   /* make sure all equations in the same node have the same status */
   /* ------------------------------------------------------------- */

   for ( inode = 0; inode < Nrows; inode+= num_PDE_eqns ) 
   {
      index3 = inode / num_PDE_eqns * num_PDE_eqns;
      status = 0;
      for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
      {
         if ( aggr_stat[index3+jnode] >= 0 ) status = 0;
      }
      if ( status == 0 ) 
      {
         procnum = aggr_stat[inode];
         for ( jnode = 1; jnode < num_PDE_eqns; jnode++ ) 
         {
            if ( aggr_stat[index3+jnode] >= 0 && procnum < 0 )
               procnum = aggr_stat[index3+jnode];
             else if ( aggr_stat[index3+jnode] >= 0 && 
                       aggr_stat[index3+jnode] < procnum )
               procnum = aggr_stat[index3+jnode];
         }
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
         {
            aggr_stat[index3+jnode] = procnum;
         }
      }
   }
        
   /* ============================================================= */
   /* send my status information to remote processors               */
   /* ============================================================= */

   for ( i = 0; i < total_send_leng; i++ ) 
   {
      int_buf[i] = aggr_stat[send_list[i]];
   }
   msgtype = 13445;
   ML_Aggregate_ExchangeData((char*)&aggr_stat[Nrows],(char*) int_buf,
      N_neighbors, neighbors,recv_leng,send_leng,msgtype,ML_INT, comm);

   /* ============================================================= */
   /* Phase 1 :                                                     */
   /*    This consists of two parts - aggregate border nodes first  */
   /*    followed by aggregating interior nodes.  This goes on      */
   /*    until all nodes are either selected or not selected.       */ 
   /* ============================================================= */

   loop_flag = 1;

   while ( loop_flag != 0 ) 
   {
      /* ========================================================== */
      /* aggregate border nodes first                               */
      /* ---------------------------------------------------------- */

      old_aggr_count = aggr_count;

      /* ---------------------------------------------------------- */
      /* look at all nodes, incremented by num_PDE_eqns             */
      /* ---------------------------------------------------------- */

      for ( inode = 0; inode < Nrows; inode+=num_PDE_eqns ) 
      {
         /* ------------------------------------------------------- */
         /* check to make sure the node is a READY node and it is   */
         /* also a border node before further considering it        */
         /* ------------------------------------------------------- */

         index3 = inode / num_PDE_eqns * num_PDE_eqns;
         count3 = 0;
         status = 1;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
         {
            if ( aggr_stat[index3+jnode] != ML_AGGR_READY ) status = 0;
            if ( ext_deg_list[index3+jnode] != 0 ) count3++;
         }
         length = 0;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
            length += int_deg_list[index3+jnode] + 
                      ext_deg_list[index3+jnode] + 1;
         length *= num_PDE_eqns;

         /* ------------------------------------------------------- */
         /* if indeed a READY and a border node, do the following   */
         /* ------------------------------------------------------- */

         if ( status == 1 && count3 > 0 )
         {
            /* ---------------------------------------------------- */
            /* first put the nodes in the supernode list            */
            /* ---------------------------------------------------- */

            supernode = (ML_SuperNode *) malloc(sizeof(ML_SuperNode));      
            supernode->list = (int*) malloc(length*sizeof(int));
            supernode->maxlength = length;
            supernode->length = num_PDE_eqns;
            for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
               supernode->list[jnode] = index3 + jnode;
            select_flag = 1;

            /* ---------------------------------------------------- */
            /* examine all of their neighbors                       */
            /* if my node is eligible to aggregate, select_flag     */
            /* will remain at 1 at the end.                         */ 
            /* ---------------------------------------------------- */

            oldleng = supernode->length;

            for ( kk = 0; kk < num_PDE_eqns; kk++ ) 
            {
               for (jnode=mat_indx[index3+kk];jnode<mat_indx[index3+kk+1];
                    jnode++) 
               {
                  index = mat_indx[jnode];
   
                  /* ---------------------------------------------- */
                  /* for each equation in my node, check to see the */
                  /* aggregation conditions are met.                */
                  /* ---------------------------------------------- */

                  index4 = index / num_PDE_eqns * num_PDE_eqns;
                  for (jj = 0; jj < num_PDE_eqns; jj++ )
                  {
                     /*-------------------------------------------- */
                     /* if my neighbors has been selected, I will   */
                     /* not use my node as seed, but rather wait    */
                     /* for others to aggregate me.                 */
                     /*-------------------------------------------- */

                     if ( aggr_stat[index4+jj] == ML_AGGR_SELECTED ) 
                     {
                        select_flag = 0;
                        procnum = mypid;

                        for (j=mat_indx[index3+kk];j<mat_indx[index3+kk+1];j++) 
                        {
                           index2 = mat_indx[j];
                           if ( index2 >= Nrows) 
                           {
                              m = mat_indx[j] - Nrows;
                              count = 0;
                              for (k = 0; k < N_neighbors; k++ ) 
                              {
                                 if ( m < (count+recv_leng[k]) ) break;
                                 count += recv_leng[k];
                              }
                              if (neighbors[k] > procnum) 
                                 procnum = neighbors[k];
                           } 
                        }

                        if (procnum==mypid) aggr_stat[index3+kk] = ML_AGGR_NOTSEL;
                        else                aggr_stat[index3+kk] = procnum;
                     } 

                     /*-------------------------------------------- */
                     /* if my neighbor is a NONSEL or READY, and it */
                     /* lives in my processor, aggregate it         */
                     /*-------------------------------------------- */

                     else if ((aggr_stat[index4+jj] == ML_AGGR_NOTSEL ||
                               aggr_stat[index4+jj] == ML_AGGR_READY) &&
                               index4+jj < Nrows)
                     {
                        if ( select_flag == 1 ) select_flag = 1;
                     }

                     /*-------------------------------------------- */
                     /* if my neighbor is a NONSEL or WAIT, and it  */
                     /* lives on another processor                  */
                     /*-------------------------------------------- */

                     else if (aggr_stat[index4+jj] >= mypid &&
                              index4+jj >= Nrows)
                     {
                        if ( select_flag == 1 ) select_flag = 1;
                     }
                     else select_flag = 0;
                  }
                  if ( select_flag == 1 ) 
                  {
                     for ( k = 0; k < num_PDE_eqns; k++ )
                     {
                        status = ML_crude_search(index4+k,supernode->length,
                                                  supernode->list);
                        if ( status < 0 )
                           supernode->list[supernode->length++] = index4+k;
                     }
                  }
               }
            }
            if ( select_flag != 1 ) supernode->length = oldleng; 

            /* ---------------------------------------------------- */
            /* synchronize all nodes in my group of num_PDE nodes   */
            /* ---------------------------------------------------- */

            procnum = aggr_stat[index3];
            for ( kk = 1; kk < num_PDE_eqns; kk++ ) 
            {
               if (aggr_stat[index3+kk] >= 0 && procnum < 0) 
                  procnum = aggr_stat[index+kk];
               else if (aggr_stat[index3+kk] >= 0 && 
                        aggr_stat[index3+kk] < procnum) 
                  procnum = aggr_stat[index+kk];
            }
            if ( procnum >= 0 )
            {
               for ( kk = 0; kk < num_PDE_eqns; kk++ ) 
                  aggr_stat[index3+kk] = procnum;
            }
            count = 0;
            for ( kk = 0; kk < num_PDE_eqns; kk++ ) 
            {
               if (aggr_stat[index3+kk] == ML_AGGR_NOTSEL ) count++;
            }
            if ( count == num_PDE_eqns )
            {
               for ( kk = 0; kk < num_PDE_eqns; kk++ ) 
                  aggr_stat[index3+kk] = ML_AGGR_NOTSEL;
            }

            /* ---------------------------------------------------- */
            /* if select_flag == 1, aggregation is successful       */
            /* ---------------------------------------------------- */

            if ( select_flag != 1 ) 
            {
               free( supernode->list );
               free( supernode );
            } 
            else 
            {
               for ( j = 0; j < supernode->length; j++ ) 
               {
                  jnode = supernode->list[j];
                  if ( jnode < exp_Nrows )
                  {
                     aggr_stat[jnode] = ML_AGGR_SELECTED;
                     aggr_index[jnode] = aggr_count;
                  }
               }
               supernode->next = NULL;
               supernode->index = aggr_count;
               if ( aggr_count == 0 ) 
               {
                  aggr_head = supernode;
                  aggr_curr = supernode;
               } else 
               {
                  aggr_curr->next = supernode;
                  aggr_curr = supernode;
               } 
               aggr_cnt_array[aggr_count++] = supernode->length;
               if ( aggr_count >= aggr_cnt_leng ) 
               {
                  itmp_array = aggr_cnt_array;
                  aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
                  nbytes = aggr_cnt_leng * sizeof( int );
                  ML_memory_alloc((void**) &aggr_cnt_array, nbytes, "AGY");
                  for ( k = 0; k < aggr_count; k++ )
                     aggr_cnt_array[k] = itmp_array[k];
                  ML_memory_free((void**) &itmp_array);
               }
            }
         }
      }
   
      /* ---------------------------------------------------------- */
      /* communicate remote node info back to remote processors     */
      /* (tell remote processor that some of their nodes have been  */
      /*  aggregated by this processor)                             */
      /* ---------------------------------------------------------- */

      msgtype = 33945 + loop_flag;
      ML_Aggregate_ExchangeData((char*)int_buf2,(char*) &(aggr_stat[Nrows]),
         N_neighbors, neighbors, send_leng, recv_leng,msgtype,ML_INT, comm);

      /* ---------------------------------------------------------- */
      /* after my processor obtains information from other          */
      /* processors about my nodes being selected, update my local  */ 
      /* node status array aggr_stat (mark ML_AGGR_SELECTED for     */
      /* local nodes that have been aggregated by remote processors)*/
      /* ---------------------------------------------------------- */

      count = 0;
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         for ( j = 0; j < send_leng[i]; j++ ) 
         {
            inode = send_list[count];
            if ( int_buf2[count] == ML_AGGR_SELECTED && 
                 aggr_stat[inode] != ML_AGGR_SELECTED ) 
            {
               aggr_stat[inode]  = ML_AGGR_SELECTED;
               aggr_index[inode] = - 100 - neighbors[i];
            } 
            count++;
         }
      }

      /* ---------------------------------------------------------- */
      /* now my aggr_stat contains latest information about the     */
      /* status of the nodes I own.  Next, send this updated info   */
      /* to other processors                                        */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < total_send_leng; i++ ) 
      {
         int_buf[i] = aggr_stat[send_list[i]];
      }
      msgtype = 13945 + loop_flag;
      ML_Aggregate_ExchangeData((char*)&aggr_stat[Nrows],(char*) int_buf,
         N_neighbors, neighbors,recv_leng,send_leng,msgtype,ML_INT, comm);

      /* ---------------------------------------------------------- */
      /* update my waiting nodes' status                            */
      /* ---------------------------------------------------------- */

      for ( inode = 0; inode < Nrows; inode++ ) 
      {
         if ( aggr_stat[inode] >= 0 )
         {
            procnum = 100000;
            for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1]; jnode++) 
            {
               index = mat_indx[jnode];
               mdiff = index - Nrows;

               if ( mdiff >= 0 )
               {
                  for ( k = 0; k <= N_neighbors; k++ )
                     if ( mdiff < sendlist_proc[k] ) break;

                  if ( aggr_stat[index] == ML_AGGR_READY )
                     procnum = neighbors[k-1];
                  else if (aggr_stat[index] >= 0 &&
                           aggr_stat[index] < mypid)
                  {
                     if ( neighbors[k-1] < procnum )
                        procnum = neighbors[k-1];
                  }
               }
            }
            if ( procnum == 100000 ) aggr_stat[inode] = ML_AGGR_READY;
            else                     aggr_stat[inode] = procnum;
         }
      }
      for ( inode = 0; inode < Nrows; inode+=num_PDE_eqns ) 
      {
         status = 1;
         index3 = inode / num_PDE_eqns * num_PDE_eqns;
         procnum = 100000;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
         {
            if ( aggr_stat[index3+jnode] >= 0 )
            {
               status = 0;
               if ( aggr_stat[index3+jnode] < procnum )
                  procnum = aggr_stat[index3+jnode];
            }
         }
         if ( status == 0 )
         {
            for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
               aggr_stat[index3+jnode] = procnum;
         }
      } 
               
      /* ---------------------------------------------------------- */
      /* now my aggr_stat contains latest information about the     */
      /* status of the nodes I own.  Next, send this updated info   */
      /* to other processors                                        */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < total_send_leng; i++ ) 
      {
         int_buf[i] = aggr_stat[send_list[i]];
      }
      msgtype = 13965 + loop_flag;
      ML_Aggregate_ExchangeData((char*)&aggr_stat[Nrows],(char*) int_buf,
         N_neighbors, neighbors,recv_leng,send_leng,msgtype,ML_INT, comm);

      /* ========================================================== */
      /* Phase 1b :                                                 */
      /*    for all internal nodes (external degree=0), see if the  */
      /*    node or its neighbors have been aggregated.  If so, go  */
      /*    to another node. If not, aggregate the node and its     */
      /*    neighbors if there are strongly coupled.                */
      /* ---------------------------------------------------------- */

      if ( ordering == 1 ) 
      {
         ML_memory_alloc((void**)&randomVector,(Nrows+1)*sizeof(int),"AGZ");
         for (i=0; i<Nrows; i++) randomVector[i] = i;
         ML_randomize(Nrows, randomVector);
      } 
      else if ( ordering == 2 )
      {
         new_node = (ML_Node *) malloc(sizeof(ML_Node));      
         new_node->node_id = 0;
         node_head = new_node;
         new_node->next = NULL;
      }
   
      for ( inode2 = 0; inode2 < Nrows; inode2++ ) 
      {
         if      ( ordering == 0 ) inode = inode2;
         else if ( ordering == 1 ) inode = randomVector[inode2];
         else if ( ordering == 2 ) 
         {
            if ( node_head == NULL ) 
            {
               for ( jnode = 0; jnode < Nrows; jnode++ ) 
               {
                  if ( aggr_stat[jnode] == ML_AGGR_READY )
                  { 
                     new_node = (ML_Node *) malloc(sizeof(ML_Node));      
                     new_node->node_id = jnode;
                     node_head = new_node;
                     new_node->next = NULL;
                     break;
                  }
               }
            }
            if ( node_head == NULL ) break;
            new_node = node_head;
            inode = new_node->node_id;
            node_head = new_node->next;
            free(new_node);
         }

         /* ------------------------------------------------------- */
         /* choose the ready nodes only (condition of Phase 1)      */
         /* ------------------------------------------------------- */

         status = 1;
         count3 = 0;
         index3 = inode / num_PDE_eqns * num_PDE_eqns;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
         {
            if ( aggr_stat[index3+jnode] == ML_AGGR_READY ) count3++;
            if ( ext_deg_list[index3+jnode] != 0 ) status = 0;
         }
         length = 0;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
         {
            length += int_deg_list[index3+jnode] + 
                      ext_deg_list[index3+jnode] + 1;
         }
         length *= num_PDE_eqns;

         /* ------------------------------------------------------- */
         /* if first condition is satified, proceed                 */
         /* ------------------------------------------------------- */

         if ( count3 == num_PDE_eqns && status == 1)
         {
            supernode = (ML_SuperNode *) malloc(sizeof(ML_SuperNode));      
            supernode->list = (int*) malloc(length*sizeof(int));
            supernode->maxlength = length;
            supernode->length = num_PDE_eqns;
            for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
               supernode->list[jnode] = index3 + jnode;
            select_flag = 1;

            /* ------------------------------------------------------- */
            /* examine all of its neighbors */
            /* ------------------------------------------------------- */

            for ( m = 0; m < num_PDE_eqns; m++ ) 
            {
               for (jnode=mat_indx[index3+m];jnode<mat_indx[index3+m+1];
                    jnode++) 
               {
                  index = mat_indx[jnode];
                  index4 = index / num_PDE_eqns * num_PDE_eqns;
                  oldleng = supernode->length;
                  for (jj = 0; jj < num_PDE_eqns; jj++ )
                  {
                     /* ---------------------------------------------- */
                     /* if my neighbor is a READY or NOTSEL, see if    */
                     /* this neighbor has already been included in the */
                     /* list. If no, include it.                       */
                     /* ---------------------------------------------- */

                     if ( aggr_stat[index4+jj] == ML_AGGR_READY || 
                          aggr_stat[index4+jj] == ML_AGGR_NOTSEL ) 
                     {
                        status = ML_crude_search(index4+jj,supernode->length,
                                                  supernode->list);
                        if ( status < 0 )
                           supernode->list[supernode->length++] = index4+jj;
                     } 

                     /* ---------------------------------------------- */
                     /* if my neighbor is SELECTED                     */
                     /* ---------------------------------------------- */

                     if ( aggr_stat[index4+jj] == ML_AGGR_SELECTED ) 
                     { 
                        select_flag = 0;
                        for (kk=0; kk<num_PDE_eqns; kk++)
                           aggr_stat[index3+kk] = ML_AGGR_NOTSEL;
                     }

                     /* ---------------------------------------------- */
                     /* if my node is not to be selected, reverse it   */
                     /* ---------------------------------------------- */

                     if (select_flag != 1) {supernode->length=oldleng;break;}
                  }
                  if ( select_flag != 1 ) break;
               }
               if ( select_flag != 1 ) break;
            }

            /* ---------------------------------------------------- */
            /* if select_flag == 1, aggregation is successful       */
            /* ---------------------------------------------------- */

            if ( select_flag != 1 ) 
            {
               free( supernode->list );
               free( supernode );
            } 
            else 
            {
               for ( j = 0; j < supernode->length; j++ ) 
               {
                  jnode = supernode->list[j];
                  if ( jnode < exp_Nrows )
                  {
                     aggr_stat[jnode] = ML_AGGR_SELECTED;
                     aggr_index[jnode] = aggr_count;
                  }
               }
               supernode->next = NULL;
               supernode->index = aggr_count;
               if ( aggr_count == 0 ) 
               {
                  aggr_head = supernode;
                  aggr_curr = supernode;
               } 
               else 
               {
                  aggr_curr->next = supernode;
                  aggr_curr = supernode;
               } 
               aggr_cnt_array[aggr_count++] = supernode->length;
               if ( aggr_count >= aggr_cnt_leng ) 
               {
                  itmp_array = aggr_cnt_array;
                  aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
                  nbytes = aggr_cnt_leng * sizeof( int );
                  ML_memory_alloc((void**) &aggr_cnt_array, nbytes, "AGa");
                  for ( k = 0; k < aggr_count; k++ )
                     aggr_cnt_array[k] = itmp_array[k];
                  ML_memory_free((void**) &itmp_array);
               }
            }
         }
      }

      /* ---------------------------------------------------------- */
      /* clean up the ordering mess                                 */
      /* ---------------------------------------------------------- */

      if ( ordering == 1 ) ML_memory_free((void**) &randomVector);
      else if ( ordering == 2 ) 
      {
         while ( node_head != NULL )
         {
            new_node = node_head;
            node_head = new_node->next;
            free( new_node );
         }
      }

      /* ---------------------------------------------------------- */
      /* send my aggr_stat information to other processors          */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < total_send_leng; i++ ) 
      {
         int_buf[i] = aggr_stat[send_list[i]];
      }
      msgtype = 38945 + loop_flag;
      ML_Aggregate_ExchangeData((char*)&aggr_stat[Nrows],(char*) int_buf,
         N_neighbors, neighbors,recv_leng,send_leng,msgtype,ML_INT, comm);

      /* ---------------------------------------------------------- */
      /* update my waiting nodes' status                            */
      /* ---------------------------------------------------------- */

      for ( inode = 0; inode < Nrows; inode++ ) 
      {
         if ( aggr_stat[inode] >= 0 )
         {
            procnum = 100000;
            for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1]; jnode++) 
            {
               index = mat_indx[jnode];
               mdiff = index - Nrows;

               if ( mdiff >= 0 )
               {
                  for ( k = 0; k <= N_neighbors; k++ )
                     if ( mdiff < sendlist_proc[k] ) break;

                  if ( aggr_stat[index] == ML_AGGR_READY ||
                       neighbors[k] < mypid )
                  {
                     if ( neighbors[k-1] < procnum )
                        procnum = aggr_stat[inode];
                  }
               }
            }
            if ( procnum == 100000 ) aggr_stat[inode] = ML_AGGR_READY;
            else                     aggr_stat[inode] = procnum;
         }
      }
      for ( inode = 0; inode < Nrows; inode+=num_PDE_eqns ) 
      {
         status = 1;
         index3 = inode / num_PDE_eqns * num_PDE_eqns;
         procnum = 100000;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
         {
            if ( aggr_stat[index3+jnode] >= 0 )
            {
               status = 0;
               if ( aggr_stat[index3+jnode] < procnum )
                  procnum = aggr_stat[index3+jnode];
            }
         }
         if ( status == 0 )
         {
            for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
               aggr_stat[index3+jnode] = procnum;
         }
      } 
               
      /* ---------------------------------------------------------- */
      /* now my aggr_stat contains latest information about the     */
      /* status of the nodes I own.  Next, send this updated info   */
      /* to other processors                                        */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < total_send_leng; i++ ) 
      {
         int_buf[i] = aggr_stat[send_list[i]];
      }
      msgtype = 23965 + loop_flag;
      ML_Aggregate_ExchangeData((char*)&aggr_stat[Nrows],(char*) int_buf,
         N_neighbors, neighbors,recv_leng,send_leng,msgtype,ML_INT, comm);

#ifdef ML_AGGR_DEBUG
      /* ---------------------------------------------------------- */
      /* output information about aggregation progress              */
      /* ---------------------------------------------------------- */

      m = 0;
      for (i = 0; i < Nrows; i++) if (aggr_stat[i] == ML_AGGR_SELECTED) m++;
      k = ML_Comm_GsumInt( comm, m);
      m = ML_Comm_GsumInt( comm, Nrows);
      j = ML_Comm_GsumInt( comm, aggr_count );
      count = 0;
      for (i = 0; i < Nrows; i++) if (aggr_stat[i] >= 0) count++;
      i = ML_Comm_GsumInt( comm, count );
      if ( mypid == 0 ) 
      {
         printf("Aggregation(C) : Phase 1 - Iteration = %d\n", loop_flag);
         printf("Aggregation(C) : Phase 1 - nodes aggregated = %d (%d)\n",k,m);
         printf("Aggregation(C) : Phase 1 - nodes waiting    = %d (%d)\n",i,m);
         printf("Aggregation(C) : Phase 1 - total aggregates = %d \n",j);
      }
#endif

      /* ---------------------------------------------------------- */
      /* check to see if further loop is needed                     */
      /* ---------------------------------------------------------- */

      count = 0;
      for ( i = 0; i < Nrows; i++ )
         if ( aggr_stat[i] >= 0 ) count++;
      count3 = ML_Comm_GsumInt( comm, count);
      if ( count3 == 0 ) loop_flag = 0;
      else
      {
         k = aggr_count - old_aggr_count;
         m = ML_Comm_GsumInt( comm, k);
         if ( m == 0 ) loop_flag = 0;
         else          loop_flag++;
      }
   }
   if ( int_buf2 != NULL ) ML_memory_free((void**) &int_buf2);
   if ( int_buf  != NULL ) ML_memory_free((void**) &int_buf);

   /* ------------------------------------------------------------- */
   /* at the end of Phase 1, all border nodes that have not been    */
   /* aggregated will be considered as NOTSEL.                      */
   /* ------------------------------------------------------------- */

   for ( inode = 0; inode < Nrows; inode++ ) 
   {
      if ( aggr_stat[inode] >= 0 ) aggr_stat[inode] = ML_AGGR_NOTSEL;
   }
   for (i = 0; i < Nrows; i++) 
      if (aggr_index[i] >= aggr_count) 
      {
         printf("WARNING (P1) : index out of range (%d,%d,%d,%d)\n",
                 mypid,i,aggr_index[i], aggr_count);
         break;
      }
#ifdef ML_MPI
   MPI_Barrier(MPI_COMM_WORLD);
#endif

   /* ============================================================= */
   /* Phase 2 :                                                     */
   /*    for all nodes, see if it can be aggregated into one of the */
   /*    existing LOCAL aggregates.                                 */
   /*    (interprocessor communication is needed here even though   */
   /*     this operation is local due to the need  to free up a     */
   /*     remote wait node)                                         */
   /* ============================================================= */

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf, nbytes, "AGb");
   else              int_buf = NULL;

   loop_flag = 1;

   while ( loop_flag != 0 ) 
   {
      count = 1;
      while ( count > 0 )
      {
         count = 0;
         for ( inode = 0; inode < Nrows; inode+=num_PDE_eqns ) 
         {
            count3 = 0;
            index3 = inode / num_PDE_eqns;
            index3 = index3 * num_PDE_eqns;
            for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
            {
               if ( aggr_stat[index3+jnode] == ML_AGGR_READY ||
                    aggr_stat[index3+jnode] == ML_AGGR_NOTSEL ) count3++;
            }

            /* ---------------------------------------------------- */
            /* if within my group all are either READY or NOTSEL    */
            /* ---------------------------------------------------- */

            search_flag = 0;
            if ( count3 == num_PDE_eqns && attach_scheme == ML_AGGR_MINRANK ) 
            {
               /* search for a neighboring aggregate that has the fewest */
               /* number of nodes                                        */

               mincount = 100000;
               for (kk = 0; kk < num_PDE_eqns; kk++ )
               {
                  for (jnode=mat_indx[index3+kk];jnode<mat_indx[index3+kk+1]; 
                       jnode++) 
                  {
                     index = mat_indx[jnode];
                     if ( index < Nrows )
                     {
                        if (aggr_stat[index] == ML_AGGR_SELECTED &&
                            aggr_index[index] >= 0) 
                        {
                           search_flag = 1;
                           m = aggr_index[index];
                           if ( aggr_cnt_array[m] < mincount ) 
                           {
                              mincount = aggr_cnt_array[m];
                              k = index;
                           }
                        }
                     }
                  }
               }
               index = aggr_index[k];
            } 
            else if ( count3 == num_PDE_eqns ) 
            {
               /* search for a neighboring aggregate that has the most   */
               /* connection to my node                                  */

               maxcount = 0;
               for (kk = 0; kk < num_PDE_eqns; kk++ )
               {
                  length = mat_indx[index3+kk+1] - mat_indx[index3+kk];
                  nbytes = length * sizeof( int );
                  if ( nbytes > 0 )
                     ML_memory_alloc((void**) &int_buf2, nbytes, "AGR");
                  length = 0;

                  for (jnode=mat_indx[index3+kk];jnode<mat_indx[index3+kk+1]; 
                       jnode++) 
                  {
                     m = mat_indx[jnode];
                     if ( aggr_index[m] >= 0 ) 
                        int_buf2[length++] = aggr_index[m];
                  }
                  if ( length > 0 ) ML_sort(length, int_buf2);
                  m = -1;
                  if ( length > 0 ) {k = int_buf2[0]; j = 1; m = k;}
                  for ( jnode = 1; jnode < length; jnode++ ) 
                  {
                     if ( int_buf2[jnode] == k ) j++;
                     else 
                     {
                        if ( j > maxcount ) 
                        {
                           maxcount = j;
                           m = k;
                        }
                        k = int_buf2[jnode];
                        j = 1;
                     }
                  }
                  if ( m >= 0 ) {search_flag = 1; index = m;}
                  if ( nbytes > 0 ) ML_memory_free((void**) &int_buf2);
               }
                  if ( mypid ==2 && index3 == 1704 ) 
                  {
                     printf("(1)HERE : index, count3 = %d %d\n", index, count3);
                     printf("(1)HERE : length = %d \n", length);
                  }
            }

            /* ------------------------------------------------- */
            /* if search is unsuccessful, see if it is a bdry    */
            /* ------------------------------------------------- */

            if ( count3 == num_PDE_eqns && search_flag == 0 ) 
            {
               bdry_flag = 0;
               for (kk = 0; kk < num_PDE_eqns; kk++ )
               {
                  for (jnode=mat_indx[index3+kk];jnode<mat_indx[index3+kk+1]; 
                       jnode++) 
                  {
                     if ((mat_indx[jnode] < index3) ||
                         (mat_indx[jnode] >= index3+num_PDE_eqns)) bdry_flag++;
                  }
               }
               if ( bdry_flag == 0 )
               {
                  nbytes = maxnnz_per_row * sizeof( int );
                  ML_memory_alloc((void**) &int_buf, nbytes, "AGA");
                  nbytes = maxnnz_per_row * sizeof( double );
                  ML_memory_alloc((void**) &dble_buf, nbytes, "AGB");
                  for (kk = 0; kk < num_PDE_eqns; kk++ )
                  {
                     nn2 = index3 + kk;
                     getrowfunc(Amatrix->data,1,&nn2,maxnnz_per_row,int_buf, 
                                    dble_buf, &m);
                     for ( j = 0; j < m; j++ ) 
                     {
                        jnode = int_buf[j];
                        if (((jnode - index3) < 0 || 
                             (jnode-index3) >= num_PDE_eqns) && 
                             jnode < Nrows && aggr_index[jnode] >= 0) 
                        {
                           index = aggr_index[jnode];
                           search_flag = 1;
                           break;
                        }
                     }
                     if ( search_flag == 1 ) break;
                  }
                  ML_memory_free((void**) &int_buf);
                  ML_memory_free((void**) &dble_buf);
               }
            }
    
            /* ------------------------------------------------- */
            /* if search is successful, put myself in it         */
            /* ------------------------------------------------- */

            if ( search_flag == 1 ) 
            {
               /*
               index = k;
               m = aggr_index[index];
               */
               m = index;
               aggr_cnt_array[m] += num_PDE_eqns;
               for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
               {
                  aggr_index[index3+jnode] = m;
                  aggr_stat[index3+jnode] = ML_AGGR_SELECTED2;
                  if ( mypid ==2 && index3 == 1704 ) 
                     printf("HERE : m, count3 = %d %d\n", m, count3);
               }
               /*###### don't need for now (too expensive)
               supernode = aggr_head;
               for ( j = 0; j < m; j++ ) supernode = supernode->next;
               if ( supernode->length >= supernode->maxlength ) 
               {
                  length = supernode->length + 1;
                  supernode->maxlength = length;
                  itmp_array = supernode->list;
                  supernode->list = (int*) malloc(length*sizeof(int));
                  for ( j = 0; j < length-1; j++ ) 
                     supernode->list[j] = itmp_array[j];
                  free( itmp_array );
               }
               supernode->list[supernode->length++] = inode;
               */
               count += num_PDE_eqns;
            } 
         }
      }
      for ( i = 0; i < exp_Nrows; i++ ) 
      {
         if (aggr_stat[i] == ML_AGGR_SELECTED2) aggr_stat[i] = ML_AGGR_SELECTED;
      }
         
      /* ---------------------------------------------------------- */
      /* communicate the information                                */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < total_send_leng; i++ ) 
      {
         int_buf[i] = aggr_stat[send_list[i]];
      }
      msgtype = 48934 + loop_flag;
      ML_Aggregate_ExchangeData((char*)&aggr_stat[Nrows],(char*) int_buf,
          N_neighbors,neighbors,recv_leng,send_leng,msgtype,ML_INT, comm);

#ifdef ML_AGGR_DEBUG
      /* ---------------------------------------------------------- */
      /* output information about aggregation progress              */
      /* ---------------------------------------------------------- */

      m = 0;
      for (i = 0; i < Nrows; i++) if (aggr_stat[i] == ML_AGGR_SELECTED) m++;
      k = ML_Comm_GsumInt( comm, m);
      m = ML_Comm_GsumInt( comm, Nrows);
      j = ML_Comm_GsumInt( comm, aggr_count );
      if ( mypid == 0 ) 
      {
         printf("Aggregation(C) : Phase 2 - Iteration = %d\n", loop_flag);
         printf("Aggregation(C) : Phase 2 - nodes aggregated = %d (%d)\n",k,m);
         printf("Aggregation(C) : Phase 2 - total aggregates = %d \n",j);
      }
#endif
      for (i = 0; i < Nrows; i++) 
         if (aggr_index[i] >= aggr_count) 
         {
            printf("WARNING (P2) : index out of range (%d,%d,%d,%d)\n",
                    mypid,i,aggr_index[i], aggr_count);
            break;
         }

      /* ---------------------------------------------------------- */
      /* checking for condition to exit the loop                    */
      /* ---------------------------------------------------------- */

      m = ML_Comm_GsumInt( comm, count );
      if ( m == 0 ) 
      { 
         loop_flag = 0;
      } else loop_flag++; 
   }
   if ( int_buf != NULL ) ML_memory_free((void**) &int_buf);

   /* ============================================================= */
   /* Phase 3 :                                                     */
   /*    for the rest of the nodes, form an aggregate for each      */
   /* ============================================================= */

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf, nbytes, "AGc");
   else              int_buf = NULL;
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf2, nbytes, "AGd");
   else              int_buf2 = NULL;

   /* ============================================================= */
   /* label all nodes that have to wait for its neighbors           */
   /* (If any of my neighbors reside in processors with processor   */
   /*  ID lower than my processor ID, then my node has to wait)     */
   /* ============================================================= */

   for ( inode = 0; inode < Nrows; inode++ ) 
   {
      if (aggr_stat[inode] != ML_AGGR_SELECTED && 
          aggr_stat[inode] != ML_AGGR_BDRY) aggr_stat[inode] = ML_AGGR_READY;
   }
   for ( inode = 0; inode < Nrows; inode++ ) 
   {
      if (ext_deg_list[inode] != 0 && aggr_stat[inode] == ML_AGGR_READY ) 
      {
         for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1]; jnode++) 
         {
            index = mat_indx[jnode];
            mdiff = index - Nrows;

            /* ---------------------------------------------------- */
            /* search for the processor the node is coming from     */
            /* ---------------------------------------------------- */

            for ( k = 0; k <= N_neighbors; k++ )
               if ( mdiff < sendlist_proc[k] ) break;

            /* ---------------------------------------------------- */
            /* if the processor number < mypid, tag it with the     */
            /* neighbor processor with the smallest rank            */
            /* ---------------------------------------------------- */

            if ( k != 0 && neighbors[k-1] < mypid && 
                 aggr_stat[index] != ML_AGGR_SELECTED )
            {
               if ( aggr_stat[inode] < 0 )
                  aggr_stat[inode] = neighbors[k-1];
               else if ( neighbors[k-1] < aggr_stat[inode] )
                  aggr_stat[inode] = neighbors[k-1];
            }
         }
      }
   }

   /* ------------------------------------------------------------- */
   /* make sure all equations in the same node have the same status */
   /* ------------------------------------------------------------- */

   if ( num_PDE_eqns > 1 )
   {
      for ( inode = 0; inode < Nrows; inode+= num_PDE_eqns ) 
      {
         index3 = inode / num_PDE_eqns * num_PDE_eqns;
         status = 0;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
         {
            if ( aggr_stat[index3+jnode] >= 0 ) status = 0;
         }
         if ( status == 0 ) 
         {
            procnum = aggr_stat[inode];
            for ( jnode = 1; jnode < num_PDE_eqns; jnode++ ) 
            {
               if ( aggr_stat[index3+jnode] >= 0 && procnum < 0 )
                  procnum = aggr_stat[index3+jnode];
                else if ( aggr_stat[index3+jnode] >= 0 && 
                          aggr_stat[index3+jnode] < procnum )
                  procnum = aggr_stat[index3+jnode];
            }
            for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
            {
               aggr_stat[index3+jnode] = procnum;
            }
         }
      }
   }

   /* ============================================================= */
   /* communicate the information                                */
   /* ---------------------------------------------------------- */

   for (i = 0; i < total_send_leng; i++) 
   {
      int_buf[i] = aggr_stat[send_list[i]];
   }
   msgtype = 48946 + loop_flag;
   ML_Aggregate_ExchangeData((char*)&aggr_stat[Nrows],(char*) int_buf,
       N_neighbors,neighbors,recv_leng,send_leng,msgtype,ML_INT, comm);

   /* ============================================================= */
   /* Phase 3 begins                                                */
   /* ------------------------------------------------------------- */

   loop_flag = 1;

   while ( loop_flag != 0 ) 
   {

      /* ========================================================== */
      /* Phase 3a : aggregate border nodes first                    */
      /* ---------------------------------------------------------- */

      old_aggr_count = aggr_count;

      for ( inode = 0; inode < Nrows; inode+=num_PDE_eqns ) 
      {
         /* ------------------------------------------------------- */
         /* check to make sure the node is a READY node and it is   */
         /* also a border node before further considering it        */
         /* ------------------------------------------------------- */

         status = 1;
         count3 = 0;
         index3 = inode / num_PDE_eqns * num_PDE_eqns;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ )
         {
            if (aggr_stat[index3+jnode] != ML_AGGR_READY) status = 0;
            if ( ext_deg_list[index3+jnode] != 0) count3++;
         } 
         length = 0;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
            length += int_deg_list[index3+jnode] + 
                      ext_deg_list[index3+jnode] + 1;
         length *= num_PDE_eqns;

         /* ------------------------------------------------------- */
         /* either a READY node, and border (group)                 */
         /* ------------------------------------------------------- */

         if ( status == 1 && count3 > 0 )
         {
            supernode = (ML_SuperNode *) malloc(sizeof(ML_SuperNode));      
            supernode->list = (int*) malloc(length*sizeof(int));
            supernode->maxlength = length;
            supernode->length = num_PDE_eqns;
            for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
               supernode->list[jnode] = index3 + jnode;
            select_flag = 1;

            /* ---------------------------------------------------- */
            /* examine all of its neighbors                         */
            /* ---------------------------------------------------- */

            for (kk = 0; kk < num_PDE_eqns; kk++ ) 
            {
               for (jnode=mat_indx[index3+kk]; jnode < mat_indx[index3+kk+1]; 
                    jnode++) 
               {
                  index = mat_indx[jnode];
                  index4 = index / num_PDE_eqns * num_PDE_eqns;
                  oldleng = supernode->length;
                  for ( jj = 0; jj < num_PDE_eqns; jj++ )
                  {
                     /*-------------------------------------------- */
                     /* if my neighbor is a READY node and it is in */
                     /* my processor, aggregate it.                 */
                     /*-------------------------------------------- */

                     if (aggr_stat[index4+jj] == ML_AGGR_READY &&
                         index4+jj < Nrows)
                     {
                        status = ML_crude_search(index4+jj,supernode->length,
                                                 supernode->list);
                        if ( status < 0 )
                           supernode->list[supernode->length++] = index4 + jj;
                     }

                     /*-------------------------------------------- */
                     /* if my neighbor is a WAIT, and it lives on   */
                     /* another processor                           */
                     /*-------------------------------------------- */

                     else if (aggr_stat[index4+jj] >= mypid &&
                              index4+jj >= Nrows)
                     {
                        status = ML_crude_search(index4+jj,supernode->length,
                                                 supernode->list);
                        if ( status < 0 )
                           supernode->list[supernode->length++] = index4 + jj;
                     }

                     /*-------------------------------------------- */
                     /* if my neighbor is a WAIT, and it lives on   */
                     /* another processor (and it's not my turn)    */
                     /*-------------------------------------------- */

                     else if (aggr_stat[index4+jj] >= 0 &&
                              aggr_stat[index4+jj] < mypid &&
                              index4+jj >= Nrows)
                     {
                        select_flag = 0;
                     }
                  }
                  if ( select_flag != 1 ) {supernode->length=oldleng;break;}
               }
            }

            /* ---------------------------------------------------- */
            /* if selected, do further processing                   */
            /* ---------------------------------------------------- */

            if ( select_flag != 1 ) 
            {
               free( supernode->list );
               free( supernode );
            } 
            else if ( select_flag == 1 && supernode->length < nullspace_dim ) 
            {
               bdry_flag = 0;
               for (kk = 0; kk < supernode->length; kk++ )
               {
                  index = supernode->list[kk];
                  nbytes = maxnnz_per_row * sizeof( int );
                  ML_memory_alloc((void**) &int_buf, nbytes, "AGA");
                  nbytes = maxnnz_per_row * sizeof( double );
                  ML_memory_alloc((void**) &dble_buf, nbytes, "AGB");
                  getrowfunc(Amatrix->data,1,&index,maxnnz_per_row,int_buf, 
                                       dble_buf, &m);
                  for ( j = 0; j < m; j++ ) 
                  {
                     jnode = int_buf[j];
                     if ( aggr_index[jnode] >= 0 ) {bdry_flag = 1; break;}
                  }
                  ML_memory_free((void**) &int_buf);
                  ML_memory_free((void**) &dble_buf);
                  if ( bdry_flag == 1 ) break;
               }
               if ( bdry_flag == 1 )
               {
                  for (kk = 0; kk < supernode->length; kk++ )
                  {
                     index = supernode->list[kk];
                     if ( index < exp_Nrows )
                     {
                        aggr_stat[index] = ML_AGGR_SELECTED;
                        aggr_index[index] = aggr_index[jnode];
                     }
                  }
                  aggr_cnt_array[aggr_index[jnode]] += supernode->length;
               }
            }
            else
            {
               for ( j = 0; j < supernode->length; j++ ) 
               {
                  jnode = supernode->list[j];
                  if ( jnode < exp_Nrows )
                  {
                     aggr_stat[jnode] = ML_AGGR_SELECTED;
                     aggr_index[jnode] = aggr_count;
                  }
               }
               supernode->next = NULL;
               supernode->index = aggr_count;
               if ( aggr_count == 0 ) 
               {
                  aggr_head = supernode;
                  aggr_curr = supernode;
               } 
               else 
               {
                  aggr_curr->next = supernode;
                  aggr_curr = supernode;
               } 
               aggr_cnt_array[aggr_count++] = supernode->length;
               if ( aggr_count >= aggr_cnt_leng ) 
               {
                  itmp_array = aggr_cnt_array;
                  aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
                  nbytes = aggr_cnt_leng * sizeof( int );
                  ML_memory_alloc((void**) &aggr_cnt_array, nbytes, "AGe");
                  for ( k = 0; k < aggr_count; k++ )
                     aggr_cnt_array[k] = itmp_array[k];
                  ML_memory_free((void**) &itmp_array);
               }
            }
         }
      }
   
      /* ---------------------------------------------------------- */
      /* communicate remote node info back to remote processors     */
      /* ---------------------------------------------------------- */

      msgtype = 34945 + loop_flag;
      ML_Aggregate_ExchangeData((char*)int_buf2,(char*)&(aggr_stat[Nrows]),
         N_neighbors,neighbors,send_leng,recv_leng, msgtype, ML_INT, comm);

      /* ---------------------------------------------------------- */
      /* combine the two information to update aggr_stat            */
      /* ---------------------------------------------------------- */

      count = 0;
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         for ( j = 0; j < send_leng[i]; j++ ) 
         {
            inode = send_list[count];
            index3 = inode / num_PDE_eqns * num_PDE_eqns;
            if ( int_buf2[count] == ML_AGGR_SELECTED && 
                 aggr_stat[inode] != ML_AGGR_SELECTED ) 
            {
               for ( jnode = 0; jnode < num_PDE_eqns; jnode++ )
               {
                  aggr_stat[index3+jnode] = ML_AGGR_SELECTED;
                  aggr_index[index3+jnode] = - 100 - neighbors[i];
               }
            } 
            count++;
         }
      }

      /* ---------------------------------------------------------- */
      /* communicate my node information to remote processors       */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < total_send_leng; i++ ) 
      {
         int_buf[i] = aggr_stat[send_list[i]];
      }
      msgtype = 14945 + loop_flag;
      ML_Aggregate_ExchangeData((char*)&aggr_stat[Nrows],(char*) int_buf,
         N_neighbors, neighbors,recv_leng,send_leng,msgtype,ML_INT,comm);

      /* ---------------------------------------------------------- */
      /* update my waiting nodes' status                            */
      /* ---------------------------------------------------------- */

      for ( inode = 0; inode < Nrows; inode++ ) 
      {
         if ( aggr_stat[inode] >= 0 )
         {
            procnum = 100000;
            for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1]; jnode++) 
            {
               index = mat_indx[jnode];
               mdiff = index - Nrows;

               if ( mdiff >= 0 )
               {
                  for ( k = 0; k <= N_neighbors; k++ )
                     if ( mdiff < sendlist_proc[k] ) break;

                  if ( aggr_stat[index] == ML_AGGR_READY )
                     procnum = neighbors[k-1];
                  else if (aggr_stat[index] >= 0 &&
                           aggr_stat[index] < mypid)
                  {
                     if ( neighbors[k-1] < procnum )
                        procnum = neighbors[k-1];
                  }
               }
            }
            if ( procnum == 100000 ) aggr_stat[inode] = ML_AGGR_READY;
            else                     aggr_stat[inode] = procnum;
         }
      }
      if ( num_PDE_eqns > 1 )
      {
         for ( inode = 0; inode < Nrows; inode+=num_PDE_eqns ) 
         {
            status = 1;
            index3 = inode / num_PDE_eqns * num_PDE_eqns;
            procnum = 100000;
            for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
            {
               if ( aggr_stat[index3+jnode] >= 0 )
               {
                  status = 0;
                  if ( aggr_stat[index3+jnode] < procnum )
                     procnum = aggr_stat[index3+jnode];
               }
            }
            if ( status == 0 )
            {
               for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
                  aggr_stat[index3+jnode] = procnum;
            }
         } 
      } 
               
      /* ---------------------------------------------------------- */
      /* now my aggr_stat contains latest information about the     */
      /* status of the nodes I own.  Next, send this updated info   */
      /* to other processors                                        */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < total_send_leng; i++ ) 
      {
         int_buf[i] = aggr_stat[send_list[i]];
      }
      msgtype = 13965 + loop_flag;
      ML_Aggregate_ExchangeData((char*)&aggr_stat[Nrows],(char*) int_buf,
         N_neighbors, neighbors,recv_leng,send_leng,msgtype,ML_INT, comm);

      /* ========================================================== */
      /* Phase 3b : aggregate interior nodes                        */
      /* ---------------------------------------------------------- */

      for ( inode = 0; inode < Nrows; inode++ ) 
      {
         /* ------------------------------------------------------- */
         /* if not considered, and is a border node                 */
         /* ------------------------------------------------------- */

         index3 = inode / num_PDE_eqns * num_PDE_eqns;
         count3 = 0;
         status = 1;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
         {
            if (aggr_stat[index3+jnode] == ML_AGGR_READY) count3++; 
            if ( ext_deg_list[index3+jnode] != 0 ) status = 0;
         }
         length = 0;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
            length += int_deg_list[index3+jnode] + 
                      ext_deg_list[index3+jnode] + 1;
         length *= num_PDE_eqns;

         /* ------------------------------------------------------- */
         /* if condition satisfied                                  */
         /* ------------------------------------------------------- */

         if ( count3 == num_PDE_eqns && status == 1 )
         {
            supernode = (ML_SuperNode *) malloc(sizeof(ML_SuperNode));      
            supernode->list = (int*) malloc(length*sizeof(int));
            supernode->maxlength = length;
            supernode->length = num_PDE_eqns;
            for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
               supernode->list[jnode] = index3 + jnode;

            /* ---------------------------------------------------- */
            /* examine all of its neighbors                         */
            /* ---------------------------------------------------- */

            for (kk = 0; kk < num_PDE_eqns; kk++) 
            {
               for (jnode = mat_indx[index3+kk]; jnode < mat_indx[index3+kk+1]; 
                    jnode++) 
               {
                  index = mat_indx[jnode];
                  index = index / num_PDE_eqns * num_PDE_eqns;
                  if ( aggr_stat[index] == ML_AGGR_READY ) 
                  {
                     status = ML_crude_search(index,supernode->length,
                                              supernode->list);
                     if ( status < 0 )
                     {
                        for (jj = 0; jj < num_PDE_eqns; jj++) 
                           supernode->list[supernode->length++] = index + jj;
                     }
                  }
               }
            }

            /* ---------------------------------------------------- */
            /* update aggregate information                         */
            /* ---------------------------------------------------- */

            for ( j = 0; j < supernode->length; j++ ) 
            {
               jnode = supernode->list[j];
               if ( jnode < exp_Nrows )
               {
                  aggr_stat[jnode] = ML_AGGR_SELECTED;
                  aggr_index[jnode] = aggr_count;
               }
            }
            supernode->next = NULL;
            supernode->index = aggr_count;
            if ( aggr_count == 0 ) 
            {
               aggr_head = supernode;
               aggr_curr = supernode;
            } 
            else 
            {
               aggr_curr->next = supernode;
               aggr_curr = supernode;
            } 
            aggr_cnt_array[aggr_count++] = supernode->length;
            if ( aggr_count >= aggr_cnt_leng ) 
            {
               itmp_array = aggr_cnt_array;
               aggr_cnt_leng = aggr_cnt_leng * 6 / 5 + 1;
               nbytes = aggr_cnt_leng * sizeof( int );
               ML_memory_alloc((void**) &aggr_cnt_array, nbytes, "AGf");
               for ( k = 0; k < aggr_count; k++ )
                  aggr_cnt_array[k] = itmp_array[k];
               ML_memory_free((void**) &itmp_array);
            }
         }
      }

      /* ---------------------------------------------------------- */
      /* send my aggr_stat information to other processors          */
      /* ---------------------------------------------------------- */

      for ( i = 0; i < total_send_leng; i++ ) 
      {
         int_buf[i] = aggr_stat[send_list[i]];
      }
      msgtype = 48945 + loop_flag;
      ML_Aggregate_ExchangeData((char*)&aggr_stat[Nrows],(char*) int_buf,
         N_neighbors, neighbors,recv_leng,send_leng,msgtype,ML_INT, comm);

      /* ---------------------------------------------------------- */
      /* update my waiting nodes' status                            */
      /* ---------------------------------------------------------- */

      for ( inode = 0; inode < Nrows; inode++ ) 
      {
         if ( aggr_stat[inode] >= 0 )
         {
            procnum = 100000;
            for (jnode=mat_indx[inode]; jnode<mat_indx[inode+1]; jnode++) 
            {
               index = mat_indx[jnode];
               mdiff = index - Nrows;

               if ( mdiff >= 0 )
               {
                  for ( k = 0; k <= N_neighbors; k++ )
                     if ( mdiff < sendlist_proc[k] ) break;

                  if ( aggr_stat[index] == ML_AGGR_READY ||
                       neighbors[k] < mypid )
                  {
                     if ( neighbors[k-1] < procnum )
                        procnum = aggr_stat[inode];
                  }
               }
            }
            if ( procnum == 100000 ) aggr_stat[inode] = ML_AGGR_READY;
            else                     aggr_stat[inode] = procnum;
         }
      }
      for ( inode = 0; inode < Nrows; inode+=num_PDE_eqns ) 
      {
         status = 1;
         index3 = inode / num_PDE_eqns * num_PDE_eqns;
         procnum = 100000;
         for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
         {
            if ( aggr_stat[index3+jnode] >= 0 )
            {
               status = 0;
               if ( aggr_stat[index3+jnode] < procnum )
                  procnum = aggr_stat[index3+jnode];
            }
         }
         if ( status == 0 )
         {
            for ( jnode = 0; jnode < num_PDE_eqns; jnode++ ) 
               aggr_stat[index3+jnode] = procnum;
         }
      } 
               
      /* ---------------------------------------------------------- */
      /* output current aggregate information                       */
      /* ---------------------------------------------------------- */

#ifdef ML_AGGR_DEBUG
      m = 0;
      for (i = 0; i < Nrows; i++) if (aggr_stat[i] == ML_AGGR_SELECTED) m++;
      k = ML_Comm_GsumInt( comm, m);
      m = ML_Comm_GsumInt( comm, Nrows);
      j = ML_Comm_GsumInt( comm, aggr_count );
      count = 0;
      for (i = 0; i < Nrows; i++) if (aggr_stat[i] == ML_AGGR_BDRY) count++;
      count = ML_Comm_GsumInt( comm, count );
      if ( mypid == 0 ) 
      {
         printf("Aggregation(C) : Phase 3 - Iteration = %d \n", loop_flag);
         printf("Aggregation(C) : Phase 3 - nodes aggregated = %d (%d)\n",k,m);
         printf("Aggregation(C) : Phase 3 - no. bdry nodes   = %d\n",count);
         printf("Aggregation(C) : Phase 3 - total aggregates = %d \n",j);
      }
#endif

      /* ---------------------------------------------------------- */
      /* check to see if further loop is needed                     */
      /* ---------------------------------------------------------- */

      count = 0;
      for ( i = 0; i < Nrows; i++ ) 
         if (aggr_stat[i]==ML_AGGR_SELECTED || aggr_stat[i]==ML_AGGR_BDRY) 
            count++;
      m = ML_Comm_GsumInt( comm, count);
      k = ML_Comm_GsumInt( comm, Nrows);
      /*k = aggr_count - old_aggr_count;
      m = ML_Comm_GsumInt( comm, k);
      */
      if ( m == k ) { loop_flag = 0;}
      else          loop_flag++;
   }
   if ( int_buf2 != NULL ) ML_memory_free((void**) &int_buf2);
   if ( int_buf  != NULL ) ML_memory_free((void**) &int_buf);
   for (i = 0; i < Nrows; i++) 
      if (aggr_index[i] >= aggr_count) 
      {
         printf("WARNING (P3) : index out of range (%d,%d,%d,%d)\n",
                 mypid,i,aggr_index[i], aggr_count);
         break;
      }

   /* ============================================================= */
   /* check for error                                               */
   /* ============================================================= */

   m = 0;
   for ( i = 0; i < Nrows; i++ ) 
      if (aggr_stat[i] == ML_AGGR_SELECTED || aggr_stat[i] == ML_AGGR_BDRY) 
         m++;
   if ( m < Nrows ) {
      printf("%d : number of unlabeled nodes = %d\n", mypid, Nrows - m);
      for ( i = 0; i < Nrows; i++ ) 
      {
         if (aggr_stat[i] != ML_AGGR_SELECTED && aggr_stat[i] != ML_AGGR_BDRY) 
            printf("%d : unlabeled node %d = %d\n", mypid, i, aggr_stat[i]);
      }
      exit(1);
#ifdef ML_AGGR_DEBUG2
      aggr_curr = aggr_head;
      for ( i = 0; i < Nrows; i++ )
         printf("   %d : node %d status = %d \n", mypid, i, aggr_stat[i]);
      while ( aggr_curr != NULL ) 
      {
         supernode = aggr_curr;
         printf("%d : aggregate number %d \n", mypid, supernode->index);
         for ( i = 0; i < supernode->length; i++ )
            printf("   %d : member = %d \n", mypid, supernode->list[i]);
         aggr_curr = aggr_curr->next;
      }   

      for ( i = 0; i < Nrows; i++ ) 
      {
         if ( aggr_stat[i] != ML_AGGR_SELECTED ) 
         {
            printf("%d : ERROR - node %d not aggregated - %d.\n", mypid, i,
                    aggr_stat[i]);
            for (jnode = mat_indx[i]; jnode < mat_indx[i+1]; jnode++) 
               printf("%d : neighbor %d(%d) = %e ===> %d %d\n",mypid,
                       mat_indx[jnode], 
                      Nrows, mat_val[jnode], aggr_index[mat_indx[jnode]],
                      aggr_stat[mat_indx[jnode]]);
         }
      }
#endif
      exit(1);
   }
  
   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */

   Ncoarse = aggr_count;

   /* ============================================================= */
   /* update aggr_index to find out which local fine node is mapped */
   /* to which coarse aggregate in remote processors                */
   /* ------------------------------------------------------------- */

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf, nbytes, "AGg");
   else              int_buf = NULL;
   nbytes = total_recv_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf2, nbytes, "AGh");
   else              int_buf2 = NULL;

   /* ------------------------------------------------------------- */
   /* send the remote node index back to remote processors, with   */
   /* added information on which remote nodes have been aggregated */
   /* by the local aggregates (and also the aggregate numbers).    */
   /* ------------------------------------------------------------- */

   offset = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      for ( j = 0; j < recv_leng[i]; j++ ) 
      {
         if ( aggr_index[Nrows+offset+j] < 0 ) int_buf2[offset+j] = -1;
         else int_buf2[offset+j] = aggr_index[Nrows+offset+j];
      }
      offset += recv_leng[i];
   }
   msgtype = 15963;
   ML_Aggregate_ExchangeData((char*) int_buf, (char*) int_buf2,
      N_neighbors, neighbors, send_leng, recv_leng, msgtype, ML_INT, comm);

   if ( int_buf2 != NULL ) ML_memory_free((void**) &int_buf2);

   /* ------------------------------------------------------------- */
   /* if int_buf[i] > 0, this means that aggr_index[send_list[i]]   */ 
   /* has been aggregated by a remote processor                     */
   /* int_buf2 is used to tabulate how many distinct aggregates     */
   /* in remote processors are used.                                */
   /* ------------------------------------------------------------- */

   offset = 0;
   m      = 0; /* store the local index offset for remote processors */ 
   new_N_recv = 0;
   nbytes = N_neighbors * sizeof(int);
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &new_recv_leng, nbytes, "AGi");
      ML_memory_alloc((void**) &new_recv_neighbors, nbytes, "AGj");
   } 
   else 
   {
      new_recv_leng = new_recv_neighbors = NULL;
   }
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      /* ---------------------------------------------------------- */
      /* find out how large an array to allocate for int_buf2       */
      /* ---------------------------------------------------------- */

      max_count = -1;
      for ( j = 0; j < send_leng[i]; j++ ) 
      {
         index = int_buf[offset+j];
         max_count = (index > max_count ) ? index : max_count;
      }
      nbytes = ( max_count + 2 ) * sizeof(int);
      if (nbytes > 0) ML_memory_alloc((void **) &int_buf2, nbytes, "AGk");

      /* ---------------------------------------------------------- */
      /* see how many distinct remote aggregates are referenced by  */
      /* local fine nodes in aggregation in proc i ==> count        */
      /* ---------------------------------------------------------- */

      for ( j = 0; j <= max_count; j++ ) int_buf2[j] = 0;
      for ( j = 0; j < send_leng[i]; j++ ) 
      {
         index = int_buf[offset+j];
         if ( index >= 0 ) int_buf2[index]++;
         if (index >= 0 && index > max_count) 
            {printf("int_buf2 error : maxcount\n");exit(1);}
      }
      count = 0;
      for ( j = 0; j <= max_count; j++ ) 
      {
         if (int_buf2[j] > 0) 
         {
            count++; int_buf2[j] = 1;
         }
      }
      for ( j = max_count; j > 0; j-- ) int_buf2[j] = int_buf2[j-1];
      int_buf2[0] = 0;
      for ( j = 0; j < max_count; j++ ) int_buf2[j+1] += int_buf2[j];

      if ( count > 0 ) 
      {
         new_recv_leng[new_N_recv] = count * nullspace_dim;
         new_recv_neighbors[new_N_recv] = neighbors[i];
         new_N_recv++;
      } 

      /* ---------------------------------------------------------- */
      /* now assign local aggregate indices to local nodes that are */
      /* aggregated by remote processors                            */
      /* ---------------------------------------------------------- */

      for ( j = 0; j < send_leng[i]; j++ ) 
      {
         index = send_list[offset+j];

         /* ------------------------------------------------------- */
         /* The first condition indicates that the local node has   */
         /* been registered to have been aggregated by remote       */
         /* aggregates.  The second condition is needed in case     */
         /* the local node is linked to more than 1 remote          */
         /* processor (but only to one aggregate though)            */
         /* int_buf2 contains local indices of remote aggregates    */
         /* ------------------------------------------------------- */

         if ( aggr_index[index] <= -100 && int_buf[offset+j] >= 0 ) 
         {
            k = int_buf[offset+j];
            aggr_index[index] = int_buf2[k] + Ncoarse + m;
         } 
      }
      if (nbytes > 0) ML_memory_free((void **) &int_buf2);
      m += count;
      offset += send_leng[i];
   }
   exp_Ncoarse = Ncoarse + m;
 
   if ( int_buf  != NULL ) ML_memory_free((void**) &int_buf);

   /* ============================================================= */
   /* check and copy aggr_index                                     */
   /* ------------------------------------------------------------- */

   level = ml_ag->cur_level;
   nbytes = Nrows * sizeof( int );
   ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), nbytes, "AGl");
   count = aggr_count;
   for ( i = 0; i < Nrows; i+=num_PDE_eqns ) 
   {
      if ( aggr_index[i] >= 0 )
      {
         for ( j = 0; j < num_PDE_eqns; j++ ) 
            ml_ag->aggr_info[level][i+j] = aggr_index[i];
         if (aggr_index[i] >= count) count = aggr_index[i] + 1;
      }
      /*else
       *{
       *   printf("%d : CoarsenCoupled error : aggr_index[%d] < 0\n",
       *          mypid,i);
       *   exit(1);
       *}*/
   }
   ml_ag->aggr_count[level] = count; /* for relaxing boundary points */
  
   /* ============================================================= */
   /* find out how many local coarse aggregates are needed by       */
   /* remote processors for interpolation (to construct the         */
   /* communicator - send info - for P)                             */
   /* ------------------------------------------------------------- */

   new_N_send = 0;
   if ( N_neighbors > 0 ) 
   {
      nbytes = N_neighbors * sizeof(int);
      ML_memory_alloc((void**) &int_buf, nbytes, "AGm");
      nbytes = Ncoarse * sizeof(int);
      ML_memory_alloc((void**) &int_buf2, nbytes, "AGn");
      for ( i = 0; i < N_neighbors; i++ ) int_buf[i] = 0;

      /* ---------------------------------------------------------- */
      /* count which remote fine nodes belong to local aggregates   */
      /* in order to generate the communication pattern for         */
      /* the interpolation operator.                                */
      /* ---------------------------------------------------------- */

      offset = Nrows; 
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         for ( j = 0; j < Ncoarse; j++ ) int_buf2[j] = 0;
         for ( j = 0; j < recv_leng[i]; j++ ) 
         {
            index = aggr_index[offset++];
            if ( index >= 0 ) int_buf2[index]++;
         }
         count = 0;
         for ( j = 0; j < Ncoarse; j++ ) if ( int_buf2[j] > 0 ) count++;
         int_buf[i] = count * nullspace_dim;
         if ( int_buf[i] > 0 ) new_N_send++;
      }

      /* ---------------------------------------------------------- */
      /* now the number of neighbors for P has been found, the next */
      /* step is to find the send_list and send_leng for the matvec */
      /* function for interpolation                                 */
      /* ---------------------------------------------------------- */

      nbytes = new_N_send * sizeof(int);
      if ( nbytes > 0 ) 
      {
         ML_memory_alloc((void**) &new_send_leng, nbytes, "AGo");
         ML_memory_alloc((void**) &new_send_neighbors, nbytes, "AGp");
         new_N_send = 0;
         for ( i = 0; i < N_neighbors; i++ ) 
         {
            if ( int_buf[i] > 0 ) 
            {
               new_send_leng[new_N_send] = int_buf[i]; 
               new_send_neighbors[new_N_send] = neighbors[i];
               new_N_send++;
            }
         }
         count = 0;
         for ( i = 0; i < new_N_send; i++ ) count += new_send_leng[i];
         nbytes = count * sizeof(int);
         ML_memory_alloc((void**) &new_send_list, nbytes, "AGq");
         offset = Nrows;
         m = count;
         count = 0;
         for ( i = 0; i < N_neighbors; i++ ) 
         {
            for ( j = 0; j < Ncoarse; j++ ) int_buf2[j] = 0;
            for ( j = 0; j < recv_leng[i]; j++ ) 
            {
               index = aggr_index[offset++];
               if ( index >= 0 ) int_buf2[index]++;
            }
            for ( j = 0; j < Ncoarse; j++ ) 
            {
               if ( int_buf2[j] > 0 ) 
               {
                  for ( jj = 0; jj < nullspace_dim; jj++ ) 
                     new_send_list[count++] = j * nullspace_dim + jj;
               } 
            } 
         } 
         if ( m != count ) 
         {
            printf("ML_Aggregate_CoarsenCoupled : internal error (1).\n");
            exit(-1);
         }
      } 
      else 
      {
         new_send_leng = NULL;
         new_send_neighbors = NULL;
         new_send_list = NULL;
      }  
      ML_memory_free((void**) &int_buf);
      ML_memory_free((void**) &int_buf2);
   } 
   else 
   {
      new_send_leng = NULL;
      new_send_neighbors = NULL;
      new_send_list = NULL;
   }

   /* ============================================================= */
   /* set up the new operator                                       */
   /* ------------------------------------------------------------- */

   new_Nrows = Nrows;
   for ( i = 0; i < new_Nrows; i++ ) 
   {
      if ( aggr_index[i] >= exp_Ncoarse ) 
      {
         printf("WARNING : index out of bound %d = %d(%d)\n", i, aggr_index[i], 
                exp_Ncoarse);
/*
         for ( j = 0; j < new_Nrows; j++ ) 
            if ( aggr_index[j] >= exp_Ncoarse )
               printf("%d : aggr_index[%5d] = %5d *\n", mypid, j, aggr_index[j]); 
            else
               printf("%d : aggr_index[%5d] = %5d \n", mypid, j, aggr_index[j]); 
         exit(1);
*/
      }
   }
   nbytes = ( new_Nrows + 1 ) * sizeof(int); 
   ML_memory_alloc((void**)&(new_ia), nbytes, "AIA");
   nbytes = new_Nrows * nullspace_dim * sizeof(int); 
   ML_memory_alloc((void**)&(new_ja), nbytes, "AJA");
   nbytes = new_Nrows * nullspace_dim * sizeof(double); 
   ML_memory_alloc((void**)&(new_val), nbytes, "AVA");
   for ( i = 0; i < new_Nrows*nullspace_dim; i++ ) new_val[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* set up the space for storing the new null space               */
   /* ------------------------------------------------------------- */

   nbytes = Ncoarse * nullspace_dim * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&(new_null),nbytes,"AGr");
   for (i = 0; i < Ncoarse*nullspace_dim*nullspace_dim; i++) 
      new_null[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* initialize the row pointer for the CSR prolongation operator  */
   /* (each row will have at most nullspace_dim nonzero entries)    */
   /* ------------------------------------------------------------- */

   for (i = 0; i <= Nrows; i++) new_ia[i] = i * nullspace_dim;

   /* ------------------------------------------------------------- */
   /* generate an array to store which aggregate has which rows.Then*/
   /* loop through the rows of A checking which aggregate each row  */
   /* is in, and adding it to the appropriate spot in rows_in_aggs  */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**)&rows_in_aggs,aggr_count*sizeof(int*),"MLs");
   for (i = 0; i < aggr_count; i++) 
   {
      rows_in_aggs[i] = (int *) malloc(aggr_cnt_array[i]*sizeof(int));
      aggr_cnt_array[i] = 0;
      if (rows_in_aggs[i] == NULL) 
      {
         printf("Error: couldn't allocate memory in CoarsenCoupled\n");
         exit(1);
      }
   }
   for (i = 0; i < exp_Nrows; i+=num_PDE_eqns) 
   {
      if ( aggr_index[i] >= 0 && aggr_index[i] < aggr_count)
      {
         for (j = 0; j < num_PDE_eqns; j++)
         {
            index = aggr_cnt_array[aggr_index[i]]++; 
            rows_in_aggs[aggr_index[i]][index] = i + j;
         }
      }
   }

   /* ------------------------------------------------------------- */
   /* allocate work arrays for QR factorization                     */
   /* work and lwork are needed for lapack's QR routine.  These     */
   /* settings seemed easiest since I don't quite understand        */
   /* what they do, but may want to do something better here later  */
   /* ------------------------------------------------------------- */

   nbytes = total_recv_leng * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&comm_val, nbytes, "AGt");
   for (i = 0; i < total_recv_leng*nullspace_dim; i++) comm_val[i] = 0.0; 
   max_agg_size = 0;
   for (i = 0; i < aggr_count; i++) 
   {
      if (aggr_cnt_array[i] > max_agg_size) max_agg_size = aggr_cnt_array[i];
   }
   nbytes = max_agg_size * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&qr_tmp, nbytes, "AGu");
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&tmp_vect, nbytes, "AGv");

   lwork  = nullspace_dim;
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&work, nbytes, "AGw");

   /* ------------------------------------------------------------- */
   /* ship the null space information to other processors           */
   /* ------------------------------------------------------------- */
 
   if (nullspace_vect != NULL) 
   {
      nbytes = total_send_leng * nullspace_dim * sizeof(double);
      ML_memory_alloc((void**) &dble_buf, nbytes,"AG1");
      nbytes = total_recv_leng * nullspace_dim * sizeof(double);
      ML_memory_alloc((void**) &dble_buf2, nbytes,"AG2");
      length = total_send_leng * nullspace_dim;
      for ( i = 0; i < total_send_leng; i++ ) 
      {
         index = send_list[i];
         for ( j = 0; j < nullspace_dim; j++ ) 
            dble_buf[i*nullspace_dim+j] = nullspace_vect[j*Nrows+index];
      }
      msgtype = 12093;
      length = sizeof(double) * nullspace_dim;
      ML_Aggregate_ExchangeData((char*)dble_buf2,(char*) dble_buf,
            N_neighbors, neighbors, recv_leng, send_leng,msgtype,length,comm);
      ML_memory_free((void**) &dble_buf);
   } 

   /* ------------------------------------------------------------- */
   /* perform block QR decomposition                                */
   /* ------------------------------------------------------------- */

   for (i = 0; i < aggr_count; i++) 
   {
      /* ---------------------------------------------------------- */
      /* set up the matrix we want to decompose into Q and R:       */
      /* ---------------------------------------------------------- */

      length = aggr_cnt_array[i];
      if (nullspace_vect == NULL) 
      {
         for (j = 0; j < length; j++)
         {
            index = rows_in_aggs[i][j];
            for (k = 0; k < nullspace_dim; k++)
            {
               if (index % num_PDE_eqns == k) qr_tmp[k*length+j] = 1.0;
               else                           qr_tmp[k*length+j] = 0.0;
            }
         }
      }
      else 
      {
         for (k = 0; k < nullspace_dim; k++)
         {
            for (j = 0; j < length; j++)
            {
               index = rows_in_aggs[i][j];
               if (index < Nrows)
               {
                  qr_tmp[k*length+j] = nullspace_vect[k*Nrows+index];
               }
               else
               {
                  qr_tmp[k*length+j] = 
                        dble_buf2[(index-Nrows)*nullspace_dim+k];
               }
            }
         }
      }

      /* ---------------------------------------------------------- */
      /* now calculate QR using an LAPACK routine                   */
      /* ---------------------------------------------------------- */

      MLFORTRAN(dgeqrf)(&(aggr_cnt_array[i]), &nullspace_dim, qr_tmp, 
                        &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
      if (info != 0)
         pr_error("Error in CoarsenCoupled : dgeqrf returned a non-zero\n");

      if (work[0] > lwork) 
      {
         lwork=(int) work[0]; 
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGx");
      }
      else lwork=work[0];
		 
      /* ---------------------------------------------------------- */
      /* the upper triangle of qr_tmp is now R, so copy that into   */
      /* the new nullspace                                          */
      /* ---------------------------------------------------------- */

      for (j = 0; j < nullspace_dim; j++)
         for (k = j; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse*nullspace_dim] = 
               qr_tmp[j+aggr_cnt_array[i]*k];
		 
      /* ---------------------------------------------------------- */
      /* to get this block of P, need to run qr_tmp through another */
      /* LAPACK function:                                           */
      /* ---------------------------------------------------------- */

      if ( aggr_cnt_array[i] < nullspace_dim )
         printf("ERROR : performing QR on a MxN matrix where M<N.\n");
      MLFORTRAN(dorgqr)(&(aggr_cnt_array[i]), &nullspace_dim, &nullspace_dim, 
              qr_tmp, &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
      if (info != 0)
         pr_error("Error in CoarsenCoupled: dorgqr returned a non-zero\n");

      if (work[0] > lwork) 
      {
         lwork=(int) work[0]; 
         ML_memory_free((void**) &work);
         ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGy");
      }
      else lwork=work[0];
		 
      /* ---------------------------------------------------------- */
      /* now copy Q over into the appropriate part of P:            */
      /* The rows of P get calculated out of order, so I assume the */
      /* Q is totally dense and use what I know of how big each Q   */
      /* will be to determine where in ia, ja, etc each nonzero in  */
      /* Q belongs.  If I did not assume this, I would have to keep */
      /* all of P in memory in order to determine where each entry  */
      /* should go                                                  */
      /* ---------------------------------------------------------- */

      for (j = 0; j < aggr_cnt_array[i]; j++)
      {
         index = rows_in_aggs[i][j];
         if ( index < Nrows )
         {
            index3 = new_ia[index];
            for (k = 0; k < nullspace_dim; k++) 
            {
               new_ja [index3+k] = i * nullspace_dim + k;
               new_val[index3+k] = qr_tmp[ k*aggr_cnt_array[i]+j];
            }
         }
         else 
         {
            index3 = (index - Nrows) * nullspace_dim;
            for (k = 0; k < nullspace_dim; k++) 
               comm_val[index3+k] = qr_tmp[ k*aggr_cnt_array[i]+j];
         }
      }
   }
	 
   ML_Aggregate_Set_NullSpace(ml_ag, num_PDE_eqns, nullspace_dim, 
                              new_null, Ncoarse*nullspace_dim);
   ML_memory_free( (void **) &new_null);
   if (nullspace_vect != NULL) ML_memory_free( (void **) &dble_buf2);

   /* ------------------------------------------------------------- */
   /* send the P rows back to its parent processor                  */
   /* ------------------------------------------------------------- */
 
   nbytes = total_send_leng * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**) &dble_buf, nbytes,"AGz");
   msgtype = 24945;
   length = sizeof(double) * nullspace_dim;
   ML_Aggregate_ExchangeData((char*)dble_buf,(char*) comm_val,
         N_neighbors, neighbors, send_leng, recv_leng,msgtype,length,comm);
   for ( i = 0; i < total_send_leng; i++ )
   {
      index = send_list[i];
      if ( aggr_index[index] >= aggr_count )
      {
         dcompare1 = 0.0;
         for ( j = 0; j < nullspace_dim; j++ )
         {
            index4 = i * nullspace_dim + j;
            dcompare1 += dble_buf[index4];
         }
         if ( dcompare1 != 0.0 )
         {
            index4 = i * nullspace_dim;
            k      = index * nullspace_dim;
            for ( j = 0; j < nullspace_dim; j++ )
            {
               new_val[k+j] = dble_buf[index4+j];
               new_ja[k+j]  = aggr_index[index]*nullspace_dim+j;
            }
         }
      }
   }
   ML_memory_free( (void **) &comm_val);
   ML_memory_free( (void **) &dble_buf);
 
   /* ------------------------------------------------------------- */
   /* check P (row sum = 1)                                         */
   /* ------------------------------------------------------------- */

   for ( i = 0; i < Nrows; i++ )
   {
      dcompare1 = 0.0;
      for (j = new_ia[i]; j < new_ia[i+1]; j++) 
      {
         dcompare1 += new_val[j];
      } 
      if ( dcompare1 == 0.0 && aggr_stat[i] != ML_AGGR_BDRY )
         printf("%d : CoarsenCoupled WARNING : rowsum(P(%d)) = 0 (%d,%d)\n",
                 mypid, i, aggr_index[i], aggr_stat[i]);
   } 
   
   /* ------------------------------------------------------------- */
   /* set up the csr_data data structure                            */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &csr_data, sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;
   (*Pmatrix) = ML_Operator_Create(comm);
   ML_Operator_Set_ApplyFuncData( *Pmatrix, nullspace_dim*Ncoarse, Nrows, 
                                  ML_EMPTY, csr_data, Nrows, NULL, 0);
   (*Pmatrix)->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   ML_memory_alloc((void**) &aggr_comm, sizeof(ML_Aggregate_Comm),"ACO");
   aggr_comm->comm = comm;
   aggr_comm->N_send_neighbors = new_N_send;
   aggr_comm->N_recv_neighbors = new_N_recv;
   aggr_comm->send_neighbors = new_send_neighbors;
   aggr_comm->recv_neighbors = new_recv_neighbors;
   aggr_comm->send_leng = new_send_leng;
   aggr_comm->recv_leng = new_recv_leng;
   aggr_comm->send_list = new_send_list;
   aggr_comm->local_nrows = Ncoarse * nullspace_dim;
   
   m = exp_Ncoarse - Ncoarse;
   ML_CommInfoOP_Generate( &((*Pmatrix)->getrow->pre_comm), 
                           ML_Aggregate_ExchangeBdry, aggr_comm, 
                           comm, Ncoarse*nullspace_dim, m*nullspace_dim);
   ML_Operator_Set_Getrow((*Pmatrix), ML_EXTERNAL, Nrows, CSR_getrows);
   ML_Operator_Set_ApplyFunc((*Pmatrix), ML_INTERNAL, CSR_matvec);
   (*Pmatrix)->max_nz_per_row = 1;

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_memory_free((void**) &comm_val);
   ML_memory_free((void**) &mat_indx);
   ML_memory_free((void**) &neighbors);
   ML_memory_free((void**) &recv_leng);
   ML_memory_free((void**) &send_leng);
   ML_memory_free((void**) &send_list);
   ML_memory_free((void**) &int_deg_list);
   ML_memory_free((void**) &ext_deg_list);
   ML_memory_free((void**) &aggr_index);
   ML_memory_free((void**) &aggr_stat);
   ML_memory_free((void**) &sendlist_proc);
   ML_memory_free((void**) &aggr_cnt_array);
   for (i = 0; i < aggr_count; i++) free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   ML_memory_free((void**)&qr_tmp);
   ML_memory_free((void**)&tmp_vect);
   ML_memory_free((void**)&work);
   if ( new_N_send > 0 ) 
   {
      ML_memory_free((void**) &new_send_leng);
      ML_memory_free((void**) &new_send_list);
      ML_memory_free((void**) &new_send_neighbors);
   }
   if ( N_neighbors > 0 ) 
   {
      ML_memory_free((void**) &new_recv_leng);
      ML_memory_free((void**) &new_recv_neighbors);
   }
   ML_memory_free((void**) &aggr_comm);
   aggr_curr = aggr_head;
   while ( aggr_curr != NULL ) 
   {
      supernode = aggr_curr;
      aggr_curr = aggr_curr->next;
      if ( supernode->length > 0 ) free( supernode->list );
      free( supernode );
   }
   return Ncoarse*nullspace_dim;
}

#ifdef ML_AGGR_PARTEST

/* ******************************************************************** */
/* Ray, what does this function do and when is it used.                 */
/*    This looks like something Dawn added. Could we leave it for now   */
/*    until I really figure out what it does?                           */
/* ******************************************************************** */

int local_to_global_row(int row)
{
   /* global_mapping is an Nrows x 2 array, where the first column
      is the global number of a row and the second column is the local
      number of that row */
				
   int i, glob_row;

   glob_row=-1;
   i=0;
   while ( i < global_nrows ) {
      if (global_mapping[i][1] == row) {
         glob_row = global_mapping[i][0];
         break;
      }
      i++;
   }
   return(glob_row);
}

#endif

/* ******************************************************************** */
/* Print information about current state of ML_Aggregate                */
/* ******************************************************************** */

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

/* ******************************************************************** */
/* Print information about operator complexity                          */
/* ******************************************************************** */

int ML_Aggregate_Print_Complexity( ML_Aggregate *ag )
{
   if ( ag->fine_complexity != 0.0 )
      printf("Smoothed Aggregation : operator complexity = %e.\n",
              ag->operator_complexity / ag->fine_complexity);
   else
      printf("Smoothed Aggregation error :  fine complexity = 0.0.\n");
   return 0;
}

/* ******************************************************************** */
/* ******************************************************************** */
/* Exchange boundary function                                           */
/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
/* ******************************************************************** */
/* Exchange data between processors given communication information     */
/* -------------------------------------------------------------------- */

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

/* ******************************************************************** */
/* destructor                                                           */
/* -------------------------------------------------------------------- */

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

int ML_Gen_Blocks_Aggregates(ML_Aggregate *ag, int level, int *nblocks, int **block_list)
{
   *nblocks = ML_Aggregate_Get_AggrCount( ag, level );
   ML_Aggregate_Get_AggrMap( ag, level, block_list);
   return 0;
}

