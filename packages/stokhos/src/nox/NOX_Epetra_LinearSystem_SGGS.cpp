
#include "NOX_Epetra_LinearSystem_SGGS.hpp"

#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "EpetraExt_BlockVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

NOX::Epetra::LinearSystemSGGS::
LinearSystemSGGS(
  Teuchos::ParameterList& printingParams, 
  Teuchos::ParameterList& linearSolverParams, 
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& det_solver_,
  const Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> >& Cijk_,
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
  const Teuchos::RCP<Epetra_Operator>& J,
  const NOX::Epetra::Vector& cloneVector,
  const NOX::Epetra::Vector& detcloneVector,
  const Teuchos::RCP<NOX::Epetra::Scaling> s):
  det_solver(det_solver_),
  Cijk(Cijk_),
  jacInterfacePtr(iJac),
  scaling(s),
  utils(printingParams)
{
  sg_op = Teuchos::rcp_dynamic_cast<Stokhos::SGOperator>(J, true);
  sg_poly = sg_op->getSGPolynomial();
  detvec = Teuchos::rcp(new Epetra_Vector(detcloneVector.getEpetraVector()));
  base_map = &(detvec->Map());
  const Epetra_BlockMap *sg_map = &(cloneVector.getEpetraVector().Map());

  MatVecTable = linearSolverParams.get("Save MatVec Table", true);
  only_use_linear = linearSolverParams.get("Only Use Linear Terms", true);
  sz = sg_poly->basis()->size();
  K_limit = sg_poly->size();
  int dim = sg_poly->basis()->dimension();
  if (only_use_linear && sg_poly->size() > dim+1)
    K_limit = dim + 1;

  sg_df_block = 
    Teuchos::rcp(new EpetraExt::BlockVector(*base_map, *sg_map));
  kx = Teuchos::rcp(new Epetra_Vector(*base_map));
  if (MatVecTable) {
    Kx_table.resize(sz);
    for (int i=0;i<sz;i++) {
      Kx_table[i].resize(K_limit);
      for (int j=0;j<K_limit;j++) {
	Kx_table[i][j] = Teuchos::rcp(new Epetra_Vector(*base_map));
      }
    }
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
  int status = sg_op->Apply(input.getEpetraVector(), 
				  result.getEpetraVector());
  return (status == 0);
}

bool NOX::Epetra::LinearSystemSGGS::
applyJacobianTranspose(const NOX::Epetra::Vector& input, 
      			      NOX::Epetra::Vector& result) const
{
  sg_op->SetUseTranspose(true);
  int status = sg_op->Apply(input.getEpetraVector(), 
			    result.getEpetraVector());
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

  const Teuchos::Array<double>& norms = sg_poly->basis()->norm_squared();
  
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

  int iter = 0;
  while (((norm_df/norm_f)>sg_tol) && (iter<max_iter)) {
    TEUCHOS_FUNC_TIME_MONITOR("Total global solve Time");
    iter++;
    
    // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
    //  ordinal_type Cijk_size = Cijk.size();
    for (int k=0; k<sz; k++) {
      df = sg_df_block->GetBlock(k);
      dx = sg_dx_block.GetBlock(k);
      df->Update(1.0, *sg_f_block.GetBlock(k), 0.0);

      int nl = Cijk->num_values(k);
      for (int l=0; l<nl; l++) {
	int i,j;
	double c;
	Cijk->value(k,l,i,j,c); 
	if (i!=0 && i<K_limit) {
	  if (MatVecTable) {
	    df->Update(-1.0*c/norms[k],*(Kx_table[j][i]),1.0);      
	  }
	  else {
	    (*sg_poly)[i].Apply(*(sg_dx_block.GetBlock(j)),*kx);
	    df->Update(-1.0*c/norms[k],*kx,1.0);      
	  }  
	}
      }
      
      (*sg_poly)[0].Apply(*dx, *kx);

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

      df->Update(-1.0, *kx, 1.0);
      
      if (MatVecTable) {
	for(int i=1;i<K_limit;i++) {
	  (*sg_poly)[i].Apply(*dx, *(Kx_table[k][i]));
	}
      }
    } //End of k loop
    
    sg_df_block->Norm2(&norm_df);
    utils.out() << "\t Gauss-Seidel relative residual norm at iteration "<< iter <<" is " << norm_df/norm_f << std::endl;
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
  EpetraExt::BlockVector sg_x_block(View, detvec->Map(), x.getEpetraVector());
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
  EpetraExt::BlockVector sg_x_block(View, detvec->Map(), x.getEpetraVector());
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
  return;
}

void NOX::Epetra::LinearSystemSGGS::
setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
}
