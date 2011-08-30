
#include "NOX_Epetra_LinearSystem_SGJacobi.hpp"

#include "Epetra_Vector.h"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "EpetraExt_BlockVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Stokhos_SGOperatorFactory.hpp"

NOX::Epetra::LinearSystemSGJacobi::
LinearSystemSGJacobi(
  Teuchos::ParameterList& printingParams, 
  Teuchos::ParameterList& linearSolverParams_, 
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& det_solver_,
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
  const Teuchos::RCP<const Stokhos::ParallelData>& sg_parallel_data,
  const Teuchos::RCP<Epetra_Operator>& J,
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<NOX::Epetra::Scaling> s):
  det_solver(det_solver_),
  epetraCijk(sg_parallel_data->getEpetraCijk()),
  jacInterfacePtr(iJac),
  base_map(base_map_),
  sg_map(sg_map_),
  scaling(s),
  utils(printingParams)
{
  sg_op = Teuchos::rcp_dynamic_cast<Stokhos::SGOperator>(J, true);
  sg_poly = sg_op->getSGPolynomial();

  sg_df_block = 
    Teuchos::rcp(new EpetraExt::BlockVector(*base_map, *sg_map));
  kx = Teuchos::rcp(new Epetra_Vector(*base_map));

  Teuchos::RCP<Teuchos::ParameterList> sgOpParams =
    Teuchos::rcp(&(linearSolverParams_.sublist("Jacobi SG Operator")), false);
  sgOpParams->set("Include Mean", false);
  if (!sgOpParams->isParameter("Only Use Linear Terms"))
    sgOpParams->set("Only Use Linear Terms", false);

  // Build new parallel Cijk if we are only using the linear terms, Cijk
  // is distributed over proc's, and Cijk includes more than just the linear
  // terms (so we have the right column map; otherwise we will be importing
  // much more than necessary)
  if (sgOpParams->get<bool>("Only Use Linear Terms") && 
      epetraCijk->isStochasticParallel()) {
    int dim = sg_basis->dimension();
    if (epetraCijk->getKEnd() > dim+1)
      epetraCijk = 
	Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(
		       *epetraCijk, 1, dim+1));
					     
  }

  Stokhos::SGOperatorFactory sg_op_factory(sgOpParams);
  Teuchos::RCP<const EpetraExt::MultiComm> sg_comm = 
    sg_parallel_data->getMultiComm();
  mat_free_op = 
    sg_op_factory.build(sg_comm, sg_basis, epetraCijk, 
			base_map, base_map, sg_map, sg_map);
}

NOX::Epetra::LinearSystemSGJacobi::
~LinearSystemSGJacobi()
{
}

bool NOX::Epetra::LinearSystemSGJacobi::
applyJacobian(const NOX::Epetra::Vector& input, 
      		     NOX::Epetra::Vector& result) const
{
  sg_op->SetUseTranspose(false);
  int status = sg_op->Apply(input.getEpetraVector(), result.getEpetraVector());

  return (status == 0);
}

bool NOX::Epetra::LinearSystemSGJacobi::
applyJacobianTranspose(const NOX::Epetra::Vector& input, 
		       NOX::Epetra::Vector& result) const
{
  sg_op->SetUseTranspose(true);
  int status = sg_op->Apply(input.getEpetraVector(), result.getEpetraVector());
  sg_op->SetUseTranspose(false);

  return (status == 0);
}

bool NOX::Epetra::LinearSystemSGJacobi::
applyJacobianInverse(Teuchos::ParameterList &params, 
		     const NOX::Epetra::Vector &input, 
		     NOX::Epetra::Vector &result)
{
  int max_iter = params.get("Max Iterations",100 );
  double sg_tol = params.get("Tolerance", 1e-12);
 
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

  Teuchos::RCP<Epetra_Vector> df, dx;
  Teuchos::ParameterList& det_solver_params = 
    params.sublist("Deterministic Solver Parameters");

  int myBlockRows = epetraCijk->numMyRows();
  int iter = 0;
  while (((norm_df/norm_f)>sg_tol) && (iter<max_iter)) {
    TEUCHOS_FUNC_TIME_MONITOR("Total global solve Time");
    iter++;

    // Compute RHS
    if (iter == 0)
      sg_df_block->Update(1.0, sg_f_block, 0.0);
    else {
      mat_free_op->Apply(sg_dx_block, *sg_df_block);
      sg_df_block->Update(1.0, sg_f_block, -1.0);
    }

    for (int i=0; i<myBlockRows; i++) {
      df = sg_df_block->GetBlock(i);
      dx = sg_dx_block.GetBlock(i);
      NOX::Epetra::Vector nox_df(df, NOX::Epetra::Vector::CreateView);
      NOX::Epetra::Vector nox_dx(dx, NOX::Epetra::Vector::CreateView);
      
      (*sg_poly)[0].Apply(*dx, *kx);
      
      dx->PutScalar(0.0);
      // Solve linear system
      {
       TEUCHOS_FUNC_TIME_MONITOR("Total deterministic solve Time");
       det_solver->applyJacobianInverse(det_solver_params, nox_df, nox_dx);
      }

      df->Update(-1.0, *kx, 1.0);
    }

    sg_df_block->Norm2(&norm_df);
    utils.out() << "\t Jacobi relative residual norm at iteration "
		<< iter <<" is " << norm_df/norm_f << std::endl;
  } //End of iter loop 

  //return status;
  return true;
}

bool NOX::Epetra::LinearSystemSGJacobi::
applyRightPreconditioning(bool useTranspose,
			  Teuchos::ParameterList& params, 
			  const NOX::Epetra::Vector& input, 
			  NOX::Epetra::Vector& result) const
{
  return false;
}

Teuchos::RCP<NOX::Epetra::Scaling> NOX::Epetra::LinearSystemSGJacobi::
getScaling()
{
  return scaling;
}

void NOX::Epetra::LinearSystemSGJacobi::
resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s)
{
  scaling = s;
  return;
}

bool NOX::Epetra::LinearSystemSGJacobi::
computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x.getEpetraVector(), 
						  *sg_op);
  sg_poly = sg_op->getSGPolynomial();
  det_solver->setJacobianOperatorForSolve(sg_poly->getCoeffPtr(0));
  mat_free_op->setupOperator(sg_poly);
  return success;
}

bool NOX::Epetra::LinearSystemSGJacobi::
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

bool NOX::Epetra::LinearSystemSGJacobi::
destroyPreconditioner() const
{
  return det_solver->destroyPreconditioner();
}

bool NOX::Epetra::LinearSystemSGJacobi::
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
NOX::Epetra::LinearSystemSGJacobi::
getPreconditionerPolicy(bool advanceReuseCounter)
{
  return det_solver->getPreconditionerPolicy(advanceReuseCounter);
} 

bool NOX::Epetra::LinearSystemSGJacobi::
isPreconditionerConstructed() const
{
  return det_solver->isPreconditionerConstructed();
}

bool NOX::Epetra::LinearSystemSGJacobi::
hasPreconditioner() const
{
  return det_solver->hasPreconditioner();
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemSGJacobi::
getJacobianOperator() const
{
  return sg_op;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemSGJacobi::
getJacobianOperator()
{
  return sg_op;
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemSGJacobi::
getGeneratedPrecOperator() const
{
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemSGJacobi::
getGeneratedPrecOperator()
{
  return Teuchos::null;
}

void NOX::Epetra::LinearSystemSGJacobi::
setJacobianOperatorForSolve(const 
      	 Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  Teuchos::RCP<const Stokhos::SGOperator> const_sg_op = 
    Teuchos::rcp_dynamic_cast<const Stokhos::SGOperator>(solveJacOp, true);
  sg_op = Teuchos::rcp_const_cast<Stokhos::SGOperator>(const_sg_op);
  sg_poly = sg_op->getSGPolynomial();
}

void NOX::Epetra::LinearSystemSGJacobi::
setPrecOperatorForSolve(const 
      	 Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
}
