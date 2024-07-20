// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Epetra_LinearSystem_MPBD.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Assert.hpp"

NOX::Epetra::LinearSystemMPBD::
LinearSystemMPBD(
  Teuchos::ParameterList& printingParams, 
  Teuchos::ParameterList& linearSolverParams, 
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& block_solver_,
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
  const Teuchos::RCP<Epetra_Operator>& J,
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<NOX::Epetra::Scaling> s):
  block_solver(block_solver_),
  jacInterfacePtr(iJac),
  base_map(base_map_),
  scaling(s),
  utils(printingParams)
{
  mp_op = Teuchos::rcp_dynamic_cast<Stokhos::BlockDiagonalOperator>(J, true);
  block_ops = mp_op->getMPOps();
  num_mp_blocks = block_ops->size();

  std::string prec_strategy = linearSolverParams.get("Preconditioner Strategy",
						     "Standard");
  if (prec_strategy == "Standard")
    precStrategy = STANDARD;
  else if (prec_strategy == "Mean")
    precStrategy = MEAN;
  else if (prec_strategy == "On the fly")
    precStrategy = ON_THE_FLY;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
		       "Invalid preconditioner strategy " << prec_strategy);

  if (precStrategy == STANDARD)
    precs.resize(num_mp_blocks);
}

NOX::Epetra::LinearSystemMPBD::
~LinearSystemMPBD()
{
}

bool NOX::Epetra::LinearSystemMPBD::
applyJacobian(const NOX::Epetra::Vector& input, 
	      NOX::Epetra::Vector& result) const
{
  mp_op->SetUseTranspose(false);
  int status = mp_op->Apply(input.getEpetraVector(), result.getEpetraVector());

  return (status == 0);
}

bool NOX::Epetra::LinearSystemMPBD::
applyJacobianTranspose(const NOX::Epetra::Vector& input, 
      			      NOX::Epetra::Vector& result) const
{
  mp_op->SetUseTranspose(true);
  int status = mp_op->Apply(input.getEpetraVector(), result.getEpetraVector());
  mp_op->SetUseTranspose(false);

  return (status == 0);
}

bool NOX::Epetra::LinearSystemMPBD::
applyJacobianInverse(Teuchos::ParameterList &params, 
		     const NOX::Epetra::Vector &input, 
		     NOX::Epetra::Vector &result)
{
  TEUCHOS_FUNC_TIME_MONITOR("Total deterministic solve Time");

  // Extract blocks
  EpetraExt::BlockVector input_block(View, *base_map, 
				     input.getEpetraVector());
  EpetraExt::BlockVector result_block(View, *base_map, 
				      result.getEpetraVector());
  result_block.PutScalar(0.0);
   
  
  Teuchos::ParameterList& block_solver_params = 
    params.sublist("Deterministic Solver Parameters");
  
  // Solve block linear systems
  bool final_status = true;
  bool status;
  for (int i=0; i<num_mp_blocks; i++) {
    NOX::Epetra::Vector nox_input(input_block.GetBlock(i), 
				  NOX::Epetra::Vector::CreateView);
    NOX::Epetra::Vector nox_result(result_block.GetBlock(i), 
				   NOX::Epetra::Vector::CreateView);
    
    block_solver->setJacobianOperatorForSolve(block_ops->getCoeffPtr(i));

    if (precStrategy == STANDARD)
      block_solver->setPrecOperatorForSolve(precs[i]);
    else if (precStrategy == ON_THE_FLY) {
      block_solver->createPreconditioner(*(prec_x->GetBlock(i)), 
					 block_solver_params, false);
    }

    status = block_solver->applyJacobianInverse(block_solver_params, nox_input, 
						nox_result);
    final_status = final_status && status;
  }

  return final_status;
}

bool NOX::Epetra::LinearSystemMPBD::
applyRightPreconditioning(bool useTranspose,
			  Teuchos::ParameterList& params, 
			  const NOX::Epetra::Vector& input, 
			  NOX::Epetra::Vector& result) const
{
  return false;
}

Teuchos::RCP<NOX::Epetra::Scaling> NOX::Epetra::LinearSystemMPBD::
getScaling()
{
  return scaling;
}

void NOX::Epetra::LinearSystemMPBD::
resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s)
{
  scaling = s;
  return;
}

bool NOX::Epetra::LinearSystemMPBD::
computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x.getEpetraVector(), 
						  *mp_op);
  block_ops = mp_op->getMPOps();
  //block_solver->setJacobianOperatorForSolve(block_ops->getCoeffPtr(0));
  return success;
}

bool NOX::Epetra::LinearSystemMPBD::
createPreconditioner(const NOX::Epetra::Vector& x, 
		     Teuchos::ParameterList& p,
		     bool recomputeGraph) const
{
  EpetraExt::BlockVector mp_x_block(View, *base_map, x.getEpetraVector());
  Teuchos::ParameterList& solverParams = 
    p.sublist("Deterministic Solver Parameters");
  bool total_success = true;
  if (precStrategy == STANDARD) {
    for (int i=0; i<num_mp_blocks; i++) {
      block_solver->setJacobianOperatorForSolve(block_ops->getCoeffPtr(i));
      if (precs[i] != Teuchos::null)
	block_solver->setPrecOperatorForSolve(precs[i]);
      bool success = 
	block_solver->createPreconditioner(*(mp_x_block.GetBlock(i)), 
					   solverParams, recomputeGraph);
      precs[i] = block_solver->getGeneratedPrecOperator();
      total_success = total_success && success;
    }
  }

  else if (precStrategy == MEAN) {
    block_solver->setJacobianOperatorForSolve(block_ops->getCoeffPtr(0));
    bool success = 
      block_solver->createPreconditioner(*(mp_x_block.GetBlock(0)), 
					 solverParams, recomputeGraph);
    total_success = total_success && success;
  }

  else if (precStrategy == ON_THE_FLY) {
    if (prec_x == Teuchos::null)
      prec_x = Teuchos::rcp(new EpetraExt::BlockVector(mp_x_block));
    else
      *prec_x = mp_x_block;
  }

  return total_success;
}

bool NOX::Epetra::LinearSystemMPBD::
destroyPreconditioner() const
{
  return block_solver->destroyPreconditioner();
}

bool NOX::Epetra::LinearSystemMPBD::
recomputePreconditioner(const NOX::Epetra::Vector& x, 
			Teuchos::ParameterList& p) const
{
  EpetraExt::BlockVector mp_x_block(View, *base_map, x.getEpetraVector());
  Teuchos::ParameterList& solverParams = 
    p.sublist("Deterministic Solver Parameters");
  bool total_success = true;

  if (precStrategy == STANDARD) {
    for (int i=0; i<num_mp_blocks; i++) {
      block_solver->setJacobianOperatorForSolve(block_ops->getCoeffPtr(i));
      if (precs[i] != Teuchos::null)
	block_solver->setPrecOperatorForSolve(precs[i]);
      bool success = 
	block_solver->recomputePreconditioner(*(mp_x_block.GetBlock(i)), 
					      solverParams);
      precs[i] = block_solver->getGeneratedPrecOperator();
      total_success = total_success && success;
    }
  }

  else if (precStrategy == MEAN) {
    block_solver->setJacobianOperatorForSolve(block_ops->getCoeffPtr(0));
    bool success = 
      block_solver->recomputePreconditioner(*(mp_x_block.GetBlock(0)), 
					    solverParams);
    total_success = total_success && success;
  }

  else if (precStrategy == ON_THE_FLY) {
    if (prec_x == Teuchos::null)
      prec_x = Teuchos::rcp(new EpetraExt::BlockVector(mp_x_block));
    else
      *prec_x = mp_x_block;
  }

  return total_success;
}

NOX::Epetra::LinearSystem::PreconditionerReusePolicyType 
NOX::Epetra::LinearSystemMPBD::
getPreconditionerPolicy(bool advanceReuseCounter)
{
  return block_solver->getPreconditionerPolicy(advanceReuseCounter);
} 

bool NOX::Epetra::LinearSystemMPBD::
isPreconditionerConstructed() const
{
  return block_solver->isPreconditionerConstructed();
}

bool NOX::Epetra::LinearSystemMPBD::
hasPreconditioner() const
{
  return block_solver->hasPreconditioner();
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemMPBD::
getJacobianOperator() const
{
  return mp_op;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemMPBD::
getJacobianOperator()
{
  return mp_op;
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemMPBD::
getGeneratedPrecOperator() const
{
  if (precStrategy == MEAN)
    return block_solver->getGeneratedPrecOperator();
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, 
      "Cannot call getGeneratedPrecOperator() unless prec strategy is Mean");
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemMPBD::
getGeneratedPrecOperator()
{
 if (precStrategy == MEAN)
    return block_solver->getGeneratedPrecOperator();
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, 
      "Cannot call getGeneratedPrecOperator() unless prec strategy is Mean");
  return Teuchos::null;
}

void NOX::Epetra::LinearSystemMPBD::
setJacobianOperatorForSolve(const 
      	 Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  Teuchos::RCP<const Stokhos::BlockDiagonalOperator> const_mp_op = 
    Teuchos::rcp_dynamic_cast<const Stokhos::BlockDiagonalOperator>(solveJacOp, 
								    true);
  mp_op = Teuchos::rcp_const_cast<Stokhos::BlockDiagonalOperator>(const_mp_op);
  block_ops = mp_op->getMPOps();
}

void NOX::Epetra::LinearSystemMPBD::
setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  if (precStrategy == MEAN)
    block_solver->setPrecOperatorForSolve(solvePrecOp);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, 
      "Cannot call setPrecOperatorForSolve() unless prec strategy is Mean");
}
