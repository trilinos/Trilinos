/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */   
/* ******************************************************************** */

/********************************************************************* */
/*          some information about the aggregates                     */
/********************************************************************* */

#ifndef __MLAGGINFO__
#define __MLAGGINFO__


#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

  extern int ML_Aggregate_Info_Setup( ML_Aggregate *ag, int MaxLevels );
  extern int ML_Aggregate_Info_Clean( ML_Aggregate *ag, int MaxLevels );
  extern int ML_Aggregate_Visualize( ML *ml, ML_Aggregate *ag, int MaxMgLevels,
				    double *x, double *y, double *z, int Ndimensions,
				    char *base_filename );
  extern int ML_Info_DomainDecomp( ML_Aggregate_Info info,
				   ML_Comm *comm, double *H, double *h );
  extern int ML_Compute_AggregateGraphRadius( int Nrows, int ia[], int ja[],
				     int dep [],
				     int *pradius, int *pNcenter );  

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef __MLAGGMETIS__ */
