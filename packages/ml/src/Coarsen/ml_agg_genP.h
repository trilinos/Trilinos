/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person    */
/* and disclaimer.                                                           */
/* ************************************************************************* */

/* ************************************************************************* */
/* data structures to hold aggregation information                           */
/* ************************************************************************* */
/* Author        : Ray Tuminaro (SNL), Charles Tong (LLNL)                   */
/* Date          : August, 1999                                              */
/* ************************************************************************* */

#ifndef __MLGENP__
#define __MLGENP__

#include "ml_common.h"
#include "ml_operator.h"
#include "ml_aggregate.h"

/* ************************************************************************* */
/* data structure to hold getrow function                                    */
/* ------------------------------------------------------------------------- */

struct ML_AGG_Matrix_Context 
{
   ML_Operator *Amat;
   double      omega;
   double      drop_tol;
   char        *near_bdry;
   int         *aggr_info;
};

#define ML_POLY_ORDER_MAX 10

struct ML_Field_Of_Values
{
  double real_max;
  double imag_max;
  double eta;
  int poly_order;
  double R_coeff[ML_POLY_ORDER_MAX];
  double P_coeff[ML_POLY_ORDER_MAX];
};


  
/* ************************************************************************* */
/* functions defined here                                                    */
/* ************************************************************************* */

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

/* ************************************************************************* */
/* functions called by users                                                 */
/* ------------------------------------------------------------------------- */

extern int ML_Gen_MGHierarchy_UsingAggregation(ML *, int start, 
                       int increment_or_decrement, ML_Aggregate *);

/* ************************************************************************* */
/* internal functions called by developers                                   */
/* ------------------------------------------------------------------------- */

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
extern int ML_AGG_Compute_Near_Bdry(ML_Operator *Amatrix, char *near_bdry);
extern int  ML_AGG_Gen_DDProlongator(ML*,int ,int,void *data,ML_Aggregate*);
extern int ML_AGG_Gen_DDProlongator2(ML *ml,int level, int clevel, void *data,
				     ML_Aggregate *ag);
extern int  ML_AGG_DD_Matvec(void *, int , double *, int, double *);
extern int  ML_AGG_DD_Getrow(void *,int, int *, int, int *, double *, int *);
extern int  ML_AGG_Extract_Diag( ML_Operator *, double *diag);
extern void ML_AGG_Matrix_Context_Clean(void *data);
extern int  ML_AGG_DD_Solve(void *data, int, double *, int, double *);
extern int  ML_AGG_Extract_Matrix(ML_Operator *mat, int *, int **, double ***);
extern int ML_AGG_Smoother_Wrapper(void *obj, int leng1, double *outvec, 
				   int leng2, double *invec);
extern int  ML_Gen_MGHierarchy_ReuseExistingOperators(ML *ml );
extern int  ML_Gen_MGHierarchy_UsingSmoothedAggr_ReuseExistingAgg(ML *ml,
								  ML_Aggregate *ag);

extern int ML_AGG_Amat_Getrows(void *data, int N_requested_rows, 
               int requested_rows[], int allocated_space, int columns[], 
               double values[], int row_lengths[]);


#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif

