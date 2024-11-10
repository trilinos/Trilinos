/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef __MLVIZSTATS__
#define __MLVIZSTATS__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*MS*/
#define ML_AGGREGATE_VIZ_STATS_ID 24680

typedef struct ML_Aggregate_Viz_Stats_Struct
{
  int id;
  double *x;
  double *y;
  double *z;
  double *material;
  int Ndim;
  int *graph_decomposition;
  int Nlocal;
  int Naggregates;
  int local_or_global;
  int is_filled;
  int MaxNodesPerAgg;
  void *Amatrix;  /* void * so that I do not have to include
		     ml_operator.h */
  /* Stuff for Zoltan */
  int zoltan_type;
  int zoltan_estimated_its;
  int zoltan_timers;
  int smoothing_steps;

} ML_Aggregate_Viz_Stats;
/*ms*/

#endif /* #ifndef __MLAGGMETIS__ */
