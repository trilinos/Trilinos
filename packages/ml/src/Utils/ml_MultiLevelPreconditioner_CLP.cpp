/*!
 *  \file ml_MultiLevelPreconditioner_CLP.cpp
 *
 *  \brief ML black-box preconditioner for Epetra_RowMatrix derived classes.
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update do Doxygen: 22-Jul-04
 *
 */


#include "ml_common.h"
#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_ifpack_wrap.h"
#include "Teuchos_ParameterList.hpp"
#include "ml_epetra_preconditioner.h"

using namespace Teuchos;
using namespace ML_Epetra;

// ================================================ ====== ==== ==== == =

#ifdef HAVE_ML_TRIUTILS
int ML_Epetra::Set(Teuchos::ParameterList & List,
		   Trilinos_Util::CommandLineParser & CLP) 
{
 
  if( CLP.Has("-ml_defaults") )
    SetDefaults(CLP.Get("-ml_defaults","DD"),List);

  // general
  if( CLP.Has("-ml_num_levels") )
    List.set("max levels",CLP.Get("-ml_num_levels",2));
  if( CLP.Has("-ml_incr_or_decr" ) )
      List.set("increasing or decreasing",CLP.Get("-ml_incr_or_decr","increasing"));
  if( CLP.Has("-ml_output" ) )
      List.set("output",CLP.Get("-ml_output",10));
  if( CLP.Has("-ml_memory" ) )
      List.set("analyze memory",CLP.Get("-ml_memory",false));
  
  // smoother
  if( CLP.Has("-ml_smoother_type") )
    List.set("smoother: type", CLP.Get("-ml_smoother_type","Gauss-Seidel"));
  
  if( CLP.Has("-ml_smoother_type_level_0") )
    List.set("smoother: type (level 0)", CLP.Get("-ml_smoother_type_level_0","Gauss-Seidel"));
  if( CLP.Has("-ml_smoother_type_level_1") )
    List.set("smoother: type (level 1)", CLP.Get("-ml_smoother_type_level_1","Gauss-Seidel"));
  if( CLP.Has("-ml_smoother_type_level_2") )
    List.set("smoother: type (level 2)", CLP.Get("-ml_smoother_type_level_2","Gauss-Seidel"));
  if( CLP.Has("-ml_smoother_type_level_3") )
    List.set("smoother: type (level 3)", CLP.Get("-ml_smoother_type_level_3","Gauss-Seidel"));
  if( CLP.Has("-ml_smoother_type_level_4") )
    List.set("smoother: type (level 4)", CLP.Get("-ml_smoother_type_level_4","Gauss-Seidel"));
  if( CLP.Has("-ml_smoother_type_level_5") )
    List.set("smoother: type (level 5)", CLP.Get("-ml_smoother_type_level_5","Gauss-Seidel"));
  if( CLP.Has("-ml_smoother_type_level_6") )
    List.set("smoother: type (level 6)", CLP.Get("-ml_smoother_type_level_6","Gauss-Seidel"));
  if( CLP.Has("-ml_smoother_type_level_7") )
    List.set("smoother: type (level 7)", CLP.Get("-ml_smoother_type_level_7","Gauss-Seidel"));
  if( CLP.Has("-ml_smoother_type_level_8") )
    List.set("smoother: type (level 8)", CLP.Get("-ml_smoother_type_level_8","Gauss-Seidel"));
  if( CLP.Has("-ml_smoother_type_level_9") )
    List.set("smoother: type (level 9)", CLP.Get("-ml_smoother_type_level_9","Gauss-Seidel"));
  
  if( CLP.Has("-ml_smoother_sweeps") )
    List.set("smoother: sweeps", CLP.Get("-ml_smoother_sweeps",1));
  if( CLP.Has("-ml_smoother_pre_or_post") )
    List.set("smoother: pre or post", CLP.Get("-ml_smoother_pre_or_post","both"));
  if( CLP.Has("-ml_smoother_damping_factor") )
    List.set("smoother: damping factor", CLP.Get("-ml_smoother_damping_factor",1.0));

  // smoother-advanced
  if( CLP.Has("-ml_RP_smoothing") )
    List.set("R and P smoothing: type", CLP.Get("-ml_RP_smoothing","classic"));
  if( CLP.Has("-ml_RP_damping") )
    List.set("R and P smoothing: damping", CLP.Get("-ml_RP_damping","fov-10"));
  if( CLP.Has("-ml_fov_not_scaled") )
    List.set("field-of-values: use diagonal scaling", false);

  
  // aggregation
  if( CLP.Has("-ml_aggr_type") )
    List.set("aggregation: type", CLP.Get("-ml_aggr_type","Uncoupled"));

  if( CLP.Has("-ml_aggr_type_level_0") )
    List.set("aggregation: type (level 0)", CLP.Get("-ml_aggr_type_level_0","Uncoupled"));
  if( CLP.Has("-ml_aggr_type_level_1") )
    List.set("aggregation: type (level 1)", CLP.Get("-ml_aggr_type_level_1","Uncoupled"));
  if( CLP.Has("-ml_aggr_type_level_2") )
    List.set("aggregation: type (level 2)", CLP.Get("-ml_aggr_type_level_2","Uncoupled"));
  if( CLP.Has("-ml_aggr_type_level_3") )
    List.set("aggregation: type (level 3)", CLP.Get("-ml_aggr_type_level_3","Uncoupled"));
  if( CLP.Has("-ml_aggr_type_level_4") )
    List.set("aggregation: type (level 4)", CLP.Get("-ml_aggr_type_level_4","Uncoupled"));
  if( CLP.Has("-ml_aggr_type_level_5") )
    List.set("aggregation: type (level 5)", CLP.Get("-ml_aggr_type_level_5","Uncoupled"));
  if( CLP.Has("-ml_aggr_type_level_6") )
    List.set("aggregation: type (level 6)", CLP.Get("-ml_aggr_type_level_6","Uncoupled"));
  if( CLP.Has("-ml_aggr_type_level_7") )
    List.set("aggregation: type (level 7)", CLP.Get("-ml_aggr_type_level_7","Uncoupled"));
  if( CLP.Has("-ml_aggr_type_level_8") )
    List.set("aggregation: type (level 8)", CLP.Get("-ml_aggr_type_level_8","Uncoupled"));
  if( CLP.Has("-ml_aggr_type_level_9") )
    List.set("aggregation: type (level 9)", CLP.Get("-ml_aggr_type_level_9","Uncoupled"));

  if( CLP.Has("-ml_num_nodes_per_aggr") )
    List.set("aggregation: nodes per aggregate", CLP.Get("-ml_num_nodes_per_aggr",512));

  if( CLP.Has("-ml_num_nodes_per_aggr_level_0") )
    List.set("aggregation: nodes per aggregate (level 0)", CLP.Get("-ml_num_nodes_per_aggr_level_0",512));
  if( CLP.Has("-ml_num_nodes_per_aggr_level_1") )
    List.set("aggregation: nodes per aggregate (level 1)", CLP.Get("-ml_num_nodes_per_aggr_level_1",512));
  if( CLP.Has("-ml_num_nodes_per_aggr_level_2") )
    List.set("aggregation: nodes per aggregate (level 2)", CLP.Get("-ml_num_nodes_per_aggr_level_2",512));
  if( CLP.Has("-ml_num_nodes_per_aggr_level_3") )
    List.set("aggregation: nodes per aggregate (level 3)", CLP.Get("-ml_num_nodes_per_aggr_level_3",512));
  if( CLP.Has("-ml_num_nodes_per_aggr_level_4") )
    List.set("aggregation: nodes per aggregate (level 4)", CLP.Get("-ml_num_nodes_per_aggr_level_4",512));
  if( CLP.Has("-ml_num_nodes_per_aggr_level_5") )
    List.set("aggregation: nodes per aggregate (level 5)", CLP.Get("-ml_num_nodes_per_aggr_level_5",512));
  if( CLP.Has("-ml_num_nodes_per_aggr_level_6") )
    List.set("aggregation: nodes per aggregate (level 6)", CLP.Get("-ml_num_nodes_per_aggr_level_6",512));
  if( CLP.Has("-ml_num_nodes_per_aggr_level_7") )
    List.set("aggregation: nodes per aggregate (level 7)", CLP.Get("-ml_num_nodes_per_aggr_level_7",512));
  if( CLP.Has("-ml_num_nodes_per_aggr_level_8") )
    List.set("aggregation: nodes per aggregate (level 8)", CLP.Get("-ml_num_nodes_per_aggr_level_8",512));
  if( CLP.Has("-ml_num_nodes_per_aggr_level_9") )
    List.set("aggregation: nodes per aggregate (level 9)", CLP.Get("-ml_num_nodes_per_aggr_level_9",512));

  if( CLP.Has("-ml_num_local_aggr") )
    List.set("aggregation: local aggregates", CLP.Get("-ml_num_local_aggr",512));

  if( CLP.Has("-ml_num_local_aggr_level_0") )
    List.set("aggregation: local aggregates (level 0)", CLP.Get("-ml_num_local_aggr_level_0",512));
  if( CLP.Has("-ml_num_local_aggr_level_1") )
    List.set("aggregation: local aggregates (level 1)", CLP.Get("-ml_num_local_aggr_level_1",512));
  if( CLP.Has("-ml_num_local_aggr_level_2") )
    List.set("aggregation: local aggregates (level 2)", CLP.Get("-ml_num_local_aggr_level_2",512));
  if( CLP.Has("-ml_num_local_aggr_level_3") )
    List.set("aggregation: local aggregates (level 3)", CLP.Get("-ml_num_local_aggr_level_3",512));
  if( CLP.Has("-ml_num_local_aggr_level_4") )
    List.set("aggregation: local aggregates (level 4)", CLP.Get("-ml_num_local_aggr_level_4",512));
  if( CLP.Has("-ml_num_local_aggr_level_5") )
    List.set("aggregation: local aggregates (level 5)", CLP.Get("-ml_num_local_aggr_level_5",512));
  if( CLP.Has("-ml_num_local_aggr_level_6") )
    List.set("aggregation: local aggregates (level 6)", CLP.Get("-ml_num_local_aggr_level_6",512));
  if( CLP.Has("-ml_num_local_aggr_level_7") )
    List.set("aggregation: local aggregates (level 7)", CLP.Get("-ml_num_local_aggr_level_7",512));
  if( CLP.Has("-ml_num_local_aggr_level_8") )
    List.set("aggregation: local aggregates (level 8)", CLP.Get("-ml_num_local_aggr_level_8",512));
  if( CLP.Has("-ml_num_local_aggr_level_9") )
    List.set("aggregation: local aggregates (level 9)", CLP.Get("-ml_num_local_aggr_level_9",512));

  if( CLP.Has("-ml_aggr_damping_factor") )
    List.set("aggregation: damping factor", CLP.Get("-ml_aggr_damping_factor",1.333));
  if( CLP.Has("-ml_compute_field_of_values") )
    List.set("aggregation: compute field of values", true);
  if( CLP.Has("-ml_compute_field_of_values_non_scaled") )
    List.set("aggregation: compute field of values for non-scaled", true);
  
  // preconditioning type
  if( CLP.Has("-ml_prec_type") ) {
    List.set("prec type", CLP.Get("-ml_prec_type","full-MGV"));
  }

  // coarse
  if( CLP.Has("-ml_coarse_type") )
    List.set("coarse: type", CLP.Get("-ml_coarse_type","Amesos-KLU"));
  if( CLP.Has("-ml_coarse_max_procs") ) 
    List.set("coarse: max processes", CLP.Get("-ml_coarse_max_procs",4));

  // eigen-analysis
  if( CLP.Has("-ml_eigen_analysis_type") )
    List.set("eigen-analysis: type", CLP.Get("-ml_eigen_analysis_type","Anorm"));
  if( CLP.Has("-ml_eigen_analysis_tol") )
    List.set("eigen-analysis: tolerance", CLP.Get("-ml_eigen_analysis_tol",1e-2));
  if( CLP.Has("-ml_compute_null_space") )
    List.set("compute null space", CLP.Has("-ml_compute_null_space"));
  if( CLP.Has("-ml_null_space_dim") )
    List.set("null space dimension", CLP.Get("-ml_null_space_dim",1));
  if( CLP.Has("-ml_add_default_null_space") )
    List.set("add default null space", CLP.Has("-ml_add_default_null_space"));

  return 0;
  
}

// ================================================ ====== ==== ==== == =

MultiLevelPreconditioner::MultiLevelPreconditioner(const Epetra_RowMatrix & RowMatrix,
						   Trilinos_Util::CommandLineParser & CLP,
						   const bool ComputePrec ) :
  RowMatrix_(&RowMatrix),
  RowMatrixAllocated_(0)
{
  Prefix_[0] = '\0';

  /* ********************************************************************** */
  /* Parse command line to get main options                                 */
  /* ********************************************************************** */

  Set(List_,CLP);
   
  /* ********************************************************************** */
  /* back to normal initialization                                          */
  /* ********************************************************************** */

  Initialize();
  
  // construct hierarchy
  if( ComputePrec == true ) ComputePreconditioner();
}
#endif

#endif /*ifdef ML_WITH_EPETRA && ML_HAVE_TEUCHOS*/

