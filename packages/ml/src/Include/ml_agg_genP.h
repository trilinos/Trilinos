/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* data structures to hold aggregation information                      */
/* ******************************************************************** */
/* Author        : Ray Tuminaro (SNL), Charles Tong (LLNL)              */
/* Organization  : Sandia National Laboratories                         */
/* Date          : August, 1999                                         */
/* ******************************************************************** */

#ifndef __MLGENP__
#define __MLGENP__

#include "ml_operator.h"
#include "ml_aggregate.h"

/* ******************************************************************** */
/* data structure to hold getrow function                               */
/* -------------------------------------------------------------------- */

struct ML_AGG_Matrix_Context 
{
   ML_Operator *Amat;
   double      omega;
   double      drop_tol;
   int         *aggr_info;
};

/* ******************************************************************** */
/* external functions                                                   */
/* -------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif
extern void ML_2matmult(ML_Operator *Mat1, ML_Operator *Mat2,
                        ML_Operator *Result);
#ifdef __cplusplus
}
#endif

/* ******************************************************************** */
/* functions defined here                                               */
/* ******************************************************************** */

#ifdef __cplusplus
extern "C" {
#endif

/* ******************************************************************** */
/* functions called by users                                            */
/* -------------------------------------------------------------------- */

extern int ML_Gen_MGHierarchy_UsingAggregation(ML *, int start, 
                       int increment_or_decrement, ML_Aggregate *);

/* ******************************************************************** */
/* internal functions called by developers                              */
/* -------------------------------------------------------------------- */

extern int ML_Gen_MGHierarchy(ML *, int fine_level,
               int (*next_level)(ML *, int, ML_Operator *, ML_Aggregate *),
               int (*user_gen_prolongator)(ML *,int,int,void *,ML_Aggregate*),
               void *data, int internal_or_external, ML_Aggregate *);
extern int  ML_AGG_Gen_Prolongator(ML*,int ,int,void *data,ML_Aggregate*);
extern int  ML_AGG_Increment_Level(ML *, int current_level, ML_Operator *Amat,
                                  ML_Aggregate *);
extern int  ML_AGG_Decrement_Level(ML *, int current_level, ML_Operator *Amat,
                                  ML_Aggregate *);
extern int  ML_AGG_Increment_Two_Level(ML *, int current_level, 
                                      ML_Operator *Amat, ML_Aggregate *);
extern int  ML_AGG_Decrement_Two_Level(ML *, int current_level, 
                                      ML_Operator *Amat, ML_Aggregate *);
extern int  ML_AGG_JacobiSmoother_Getrows(void *data, int N_requested_rows,
               int requested_rows[], int allocated_space, int columns[],
               double values[], int row_lengths[]);
extern int  ML_AGG_Gen_DDProlongator(ML*,int ,int,void *data,ML_Aggregate*);
extern int  ML_AGG_DD_Matvec(void *, int , double *, int, double *);
extern int  ML_AGG_DD_Getrow(void *,int, int *, int, int *, double *, int *);
extern int  ML_AGG_Extract_Diag( ML_Operator *, double *diag);
extern void ML_AGG_Matrix_Context_Clean(void *data);
extern int  ML_AGG_DD_Solve(void *data, int, double *, int, double *);
extern int  ML_AGG_Extract_Matrix(ML_Operator *mat, int *, int **, double ***);

#ifdef __cplusplus
}
#endif

#endif

