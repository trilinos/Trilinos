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

#include "ml_defs.h"
#include "ml_struct.h"
#include "ml_agg_genP.h"
#include "ml_aggregate.h"
#include "ml_bicgstabl.h"
#include "ml_cg.h"
#include "ml_comm.h"
#include "ml_gmres.h"
#include "ml_grid.h"
#include "ml_gridfunc.h"
#include "ml_krylov.h"
#include "ml_operator.h"
#include "mli_solver.h"
#include "ml_pde.h"
#include "ml_solver.h"
#include "ml_vec.h"

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

