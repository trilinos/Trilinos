/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Include file                                                         */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : January, 1998                                        */
/* ******************************************************************** */

#ifndef _MLINCLUDE_
#define _MLINCLUDE_

#include "Include/ml_defs.h"
#include "Main/ml_struct.h"
#include "Comm/ml_comm.h"
#include "Utils/ml_memory.h"
#include "Operator/ml_mat_formats.h"
#include "Operator/ml_rap.h"
#include "Utils/ml_vec.h"
#include "Utils/ml_intlist.h"  
#include "FEGrid/ml_elementagx.h"
#include "Smoother/ml_solver.h"
#include "FEGrid/ml_gridfunc.h"
#include "Operator/ml_operatoragx.h"
#include "FEGrid/ml_gridagx.h"
#include "Comm/ml_comminfoagx.h"
#include "Utils/ml_utils.h"
#include "Utils/ml_aztec_utils.h"
#include "Coarsen/ml_ggraph.h"
#include "Coarsen/ml_aggregate.h"
#include "Operator/ml_op_utils.h"
#include "Coarsen/ml_agg_genP.h"
#include "Krylov/ml_krylov.h"
#include "Krylov/ml_cg.h"
#include "Krylov/ml_gmres.h"
#include "Krylov/ml_bicgstabl.h"
#include "Main/mli_solver.h"
#include "Utils/ml_rbm.h"

#ifdef __cplusplus
extern "C"
{
#endif

extern int  ML_compute_basis_coefficients2D(void *grid, int element_index, 
               double *set_of_coord,int ncoord,double *coefs,int *coef_ptr);

extern int  ML_compute_basis_coefficients3D(void *grid, int element_index, 
               double *set_of_coord,int ncoord,double *coefs,int *coef_ptr);

/*
extern void ML_setup_grid_xsfer_op(void *f_grid, ML_GridFunc *fgrid_fcns,
               void *c_grid, ML_GridFunc *cgrid_fcns, void **xsfer, 
               ML_Comm *comm );
*/

#ifdef __cplusplus
}
#endif

#endif

