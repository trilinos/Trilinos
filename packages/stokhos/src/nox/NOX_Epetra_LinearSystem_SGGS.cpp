// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Epetra_LinearSystem_SGGS.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockUtility.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

NOX::Epetra::LinearSystemSGGS::
LinearSystemSGGS(
  Teuchos::ParameterList& printingParams, 
  Teuchos::ParameterList& linearSolverParams, 
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& det_solver_,
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::ParallelData>& sg_parallel_data,
  const Teuchos::RCP<Epetra_Operator>& J,
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<NOX::Epetra::Scaling> s):
  det_solver(det_solver_),
  sg_basis(sg_basis_),
  epetraCijk(sg_parallel_data->getEpetraCijk()),
  is_stoch_parallel(epetraCijk->isStochasticParallel()),
  stoch_row_map(epetraCijk->getStochasticRowMap()),
  Cijk(epetraCijk->getParallelCijk()),
  jacInterfacePtr(iJac),
  base_map(base_map_),
  sg_map(sg_map_),
  scaling(s),
  utils(printingParams),
  is_parallel(epetraCijk->isStochasticParallel())
{
  sg_op = Teuchos::rcp_dynamic_cast<Stokhos::SGOperator>(J, true);
  sg_poly = sg_op->getSGPolynomial();

  sg_df_block = Teuchos::rcp(new EpetraExt::BlockVector(*base_map, *sg_map));
  sg_y_block = Teuchos::rcp(new EpetraExt::BlockVector(*base_map, *sg_map));
  kx = Teuchos::rcp(new Epetra_Vector(*base_map));

  only_use_linear = linearSolverParams.get("Only Use Linear Terms", false);
  k_limit = sg_poly->size();
  int dim = sg_basis->dimension();
  if (only_use_linear && sg_poly->size() > dim+1)
    k_limit = dim + 1;

  if (is_parallel) {
    Teuchos::RCP<const Epetra_BlockMap> stoch_col_map =
      epetraCijk->getStochasticColMap();
    sg_col_map =
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(*base_map,
							     *stoch_col_map,
							     sg_map->Comm()));
    col_exporter = Teuchos::rcp(new Epetra_Export(*sg_col_map, *sg_map));
    sg_df_col = Teuchos::rcp(new EpetraExt::BlockVector(*base_map, 
							*sg_col_map));
    sg_df_tmp = Teuchos::rcp(new EpetraExt::BlockVector(*base_map, *sg_map));
  }
}

NOX::Epetra::LinearSystemSGGS::
~LinearSystemSGGS()
{
}

bool NOX::Epetra::LinearSystemSGGS::
applyJacobian(const NOX::Epetra::Vector& input, 
	      NOX::Epetra::Vector& result) const
{
  sg_op->SetUseTranspose(false);
  int status = sg_op->Apply(input.getEpetraVector(), result.getEpetraVector());

  return (status == 0);
}

bool NOX::Epetra::LinearSystemSGGS::
applyJacobianTranspose(const NOX::Epetra::Vector& input, 
      			      NOX::Epetra::Vector& result) const
{
  sg_op->SetUseTranspose(true);
  int status = sg_op->Apply(input.getEpetraVector(), result.getEpetraVector());
  sg_op->SetUseTranspose(false);

  return (status == 0);
}

bool NOX::Epetra::LinearSystemSGGS::
applyJacobianInverse(Teuchos::ParameterList &params, 
		     const NOX::Epetra::Vector &input, 
		     NOX::Epetra::Vector &result)
{
  int max_iter = params.get("Max Iterations",100 );
  double sg_tol = params.get("Tolerance", 1e-12);
  bool scale_op = params.get("Scale Operator by Inverse Basis Norms", true);
  
  // Extract blocks
  EpetraExt::BlockVector sg_dx_block(View, *base_map, result.getEpetraVector());
  EpetraExt::BlockVector sg_f_block(View, *base_map, input.getEpetraVector());
  sg_dx_block.PutScalar(0.0);
   
  // Compute initial residual norm
  double norm_f,norm_df;
  sg_f_block.Norm2(&norm_f);
  sg_op->Apply(sg_dx_block, *sg_df_block);
  sg_df_block->Update(-1.0, sg_f_block, 1.0);
  sg_df_block->Norm2(&norm_df);
  
  Teuchos::RCP<Epetra_Vector> f, df, dx;
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();

  // sg_df_block stores the right-hand-sides for the Gauss-Seidel iteration
  // sg_y_block stores the residual, i.e., A*x-b
  // sg_df_col stores off-processor contributions to the RHS and residual
  // In parallel, this is essentially domain decomposition using Gauss-Seidel
  // in each domain and Jacobi across domains.

  int iter = 0;
  while (((norm_df/norm_f)>sg_tol) && (iter<max_iter)) {
    TEUCHOS_FUNC_TIME_MONITOR("Total global solve Time");
    iter++;
    
    sg_y_block->Update(1.0, sg_f_block, 0.0);
    if (is_parallel)
      sg_df_col->PutScalar(0.0);
    
    for (Cijk_type::i_iterator i_it=Cijk->i_begin(); 
	 i_it!=Cijk->i_end(); ++i_it) {
      int i = Stokhos::index(i_it);
      f = sg_f_block.GetBlock(i);
      df = sg_df_block->GetBlock(i);
      dx = sg_dx_block.GetBlock(i);

      dx->PutScalar(0.0);
      Teuchos::ParameterList& det_solver_params = 
	params.sublist("Deterministic Solver Parameters");
      NOX::Epetra::Vector nox_df(df, NOX::Epetra::Vector::CreateView);
      NOX::Epetra::Vector nox_dx(dx, NOX::Epetra::Vector::CreateView);
      // Solve linear system
      {
	TEUCHOS_FUNC_TIME_MONITOR("Total deterministic solve Time");
	det_solver->applyJacobianInverse(det_solver_params, nox_df, nox_dx);
      }

      df->Update(1.0, *f, 0.0);
      
      for (Cijk_type::ik_iterator k_it = Cijk->k_begin(i_it);
	   k_it != Cijk->k_end(i_it); ++k_it) {
	int k = Stokhos::index(k_it);
	if (k!=0 && k<k_limit) {
	  (*sg_poly)[k].Apply(*dx, *kx);
	  for (Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	       j_it != Cijk->j_end(k_it); ++j_it) {
	    int j = Stokhos::index(j_it);
	    int j_gid = epetraCijk->GCID(j);
	    double c = Stokhos::value(j_it);
	    if (scale_op)
	      c /= norms[j_gid];
	    bool owned = epetraCijk->myGRID(j_gid);
	    if (!is_parallel || owned) {
	      sg_df_block->GetBlock(j)->Update(-c, *kx, 1.0);
	      sg_y_block->GetBlock(j)->Update(-c, *kx, 1.0);
	    }
	    else
	      sg_df_col->GetBlock(j)->Update(-c, *kx, 1.0);
	  }
	}
      }
      
      (*sg_poly)[0].Apply(*dx, *kx);
      sg_y_block->GetBlock(i)->Update(-1.0, *kx, 1.0);
      
    } //End of k loop

    // Add in contributions from off-process
    if (is_parallel) {
      sg_df_tmp->Export(*sg_df_col, *col_exporter, InsertAdd);
      sg_df_block->Update(1.0, *sg_df_tmp, 1.0);
      sg_y_block->Update(1.0, *sg_df_tmp, 1.0);
    }

    sg_y_block->Norm2(&norm_df);
    utils.out() << "\t Gauss-Seidel relative residual norm at iteration "
		<< iter <<" is " << norm_df/norm_f << std::endl;
  } //End of iter loop

  //return status;
  return true;
}

bool NOX::Epetra::LinearSystemSGGS::
applyRightPreconditioning(bool useTranspose,
			  Teuchos::ParameterList& params, 
			  const NOX::Epetra::Vector& input, 
			  NOX::Epetra::Vector& result) const
{
  return false;
}

Teuchos::RCP<NOX::Epetra::Scaling> NOX::Epetra::LinearSystemSGGS::
getScaling()
{
  return scaling;
}

void NOX::Epetra::LinearSystemSGGS::
resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s)
{
  scaling = s;
  return;
}

bool NOX::Epetra::LinearSystemSGGS::
computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x.getEpetraVector(), 
						  *sg_op);
  sg_poly = sg_op->getSGPolynomial();
  det_solver->setJacobianOperatorForSolve(sg_poly->getCoeffPtr(0));
  return success;
}

bool NOX::Epetra::LinearSystemSGGS::
createPreconditioner(const NOX::Epetra::Vector& x, 
		     Teuchos::ParameterList& p,
		     bool recomputeGraph) const
{
  EpetraExt::BlockVector sg_x_block(View, *base_map, x.getEpetraVector());
  bool success = 
    det_solver->createPreconditioner(*(sg_x_block.GetBlock(0)), 
				   p.sublist("Deterministic Solver Parameters"),
				   recomputeGraph);

  return success;
}

bool NOX::Epetra::LinearSystemSGGS::
destroyPreconditioner() const
{
  return det_solver->destroyPreconditioner();
}

bool NOX::Epetra::LinearSystemSGGS::
recomputePreconditioner(const NOX::Epetra::Vector& x, 
      		Teuchos::ParameterList& linearSolverParams) const
{  
  EpetraExt::BlockVector sg_x_block(View, *base_map, x.getEpetraVector());
  bool success = 
    det_solver->recomputePreconditioner(
      *(sg_x_block.GetBlock(0)), 
      linearSolverParams.sublist("Deterministic Solver Parameters"));

  return success;

}

NOX::Epetra::LinearSystem::PreconditionerReusePolicyType 
NOX::Epetra::LinearSystemSGGS::
getPreconditionerPolicy(bool advanceReuseCounter)
{
  return det_solver->getPreconditionerPolicy(advanceReuseCounter);
} 

bool NOX::Epetra::LinearSystemSGGS::
isPreconditionerConstructed() const
{
  return det_solver->isPreconditionerConstructed();
}

bool NOX::Epetra::LinearSystemSGGS::
hasPreconditioner() const
{
  return det_solver->hasPreconditioner();
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemSGGS::
getJacobianOperator() const
{
  return sg_op;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemSGGS::
getJacobianOperator()
{
  return sg_op;
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemSGGS::
getGeneratedPrecOperator() const
{
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemSGGS::
getGeneratedPrecOperator()
{
  return Teuchos::null;
}

void NOX::Epetra::LinearSystemSGGS::
setJacobianOperatorForSolve(const 
      	 Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  Teuchos::RCP<const Stokhos::SGOperator> const_sg_op = 
    Teuchos::rcp_dynamic_cast<const Stokhos::SGOperator>(solveJacOp, true);
  sg_op = Teuchos::rcp_const_cast<Stokhos::SGOperator>(const_sg_op);
  sg_poly = sg_op->getSGPolynomial();
}

void NOX::Epetra::LinearSystemSGGS::
setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
}
