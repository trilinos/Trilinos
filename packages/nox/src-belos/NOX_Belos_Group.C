// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_Belos_Group.H"
#include "NOX_Belos_MultiVector.H"
#include "NOX_Belos_JacobianOperator.H"
#include "NOX_Belos_PreconditionOperator.H"
#include "NOX_Parameter_List.H"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblemManager.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestMaxRestarts.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosBlockGmres.hpp"
#include "BelosBlockCG.hpp"

NOX::Belos::Group::Group(NOX::Abstract::Group& g,
			 NOX::Parameter::List& printParams)
  : grpPtr(&g),
    ownsGroup(false),
    newtonVecPtr(g.getX().clone(NOX::ShapeCopy)),
    isValidNewton(false),
    utils(printParams),
    myPID(printParams.getParameter("MyPID", 0))
{
}

NOX::Belos::Group::Group(const NOX::Belos::Group& source, 
			  NOX::CopyType type)
  : grpPtr(source.grpPtr->clone(type)),
    ownsGroup(true),
    newtonVecPtr(source.newtonVecPtr->clone(type)),
    isValidNewton(false),
    utils(source.utils),
    myPID(source.myPID)
{
  if (type == NOX::DeepCopy)
    isValidNewton = source.isValidNewton;
}


NOX::Belos::Group::~Group() 
{
  if (ownsGroup)
    delete grpPtr;
  delete newtonVecPtr;
}

NOX::Belos::Group&
NOX::Belos::Group::operator=(const NOX::Belos::Group& source) 
{

  // Protect against A = A
  if (this != &source) {
    
    // Copy values
    *grpPtr = *source.grpPtr;
    *newtonVecPtr = *source.newtonVecPtr;
    isValidNewton = source.isValidNewton;
    utils = source.utils;
    myPID = source.myPID;

  }

  return *this;
}

NOX::Abstract::Group&
NOX::Belos::Group::operator=(const NOX::Abstract::Group& source)
{
  return *this = 
    dynamic_cast<const NOX::Belos::Group&>(source);
}

NOX::Abstract::Group*
NOX::Belos::Group::clone(NOX::CopyType type) const 
{
  return new NOX::Belos::Group(*this, type);
}

void
NOX::Belos::Group::setX(const NOX::Abstract::Vector& y) 
{
  grpPtr->setX(y);
  resetIsValid();
}

void
NOX::Belos::Group::computeX(const NOX::Abstract::Group& g, 
			     const NOX::Abstract::Vector& d,
			     double step) 
{
  const NOX::Belos::Group& belos_g = dynamic_cast<const NOX::Belos::Group&>(g);
  grpPtr->computeX(*(belos_g.grpPtr),d,step);
  resetIsValid();
}

NOX::Abstract::Group::ReturnType
NOX::Belos::Group::computeF() 
{
  return grpPtr->computeF();
}

NOX::Abstract::Group::ReturnType
NOX::Belos::Group::computeJacobian() 
{
  return grpPtr->computeJacobian();
}

NOX::Abstract::Group::ReturnType
NOX::Belos::Group::computeGradient() 
{
  return grpPtr->computeGradient();
}
   
NOX::Abstract::Group::ReturnType
NOX::Belos::Group::computeNewton(NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  if (!isF()) {
    cerr << "ERROR: NOX::Belos::Group::computeNewton() - invalid RHS" << endl;
    throw "NOX Error";
  }

  if (!isJacobian()) {
    cerr << "ERROR: NOX::Belos::Group::computeNewton() - invalid Jacobian" 
	 << endl;
    throw "NOX Error";
  }

  Abstract::Group::ReturnType status;

  // zero out newton vec -- used as initial guess for some linear solvers
  newtonVecPtr->init(0.0);

  status = applyJacobianInverse(params, getF(), *newtonVecPtr);
 
  newtonVecPtr->scale(-1.0);

  // Update state EVEN IF LINEAR SOLVE FAILED
  // We still may want to use the vector even it it just missed it's 
  isValidNewton = true;

  return status;
}

NOX::Abstract::Group::ReturnType
NOX::Belos::Group::applyJacobian(const NOX::Abstract::Vector& input,
				  NOX::Abstract::Vector& result) const 
{
  return grpPtr->applyJacobian(input, result);
}

NOX::Abstract::Group::ReturnType
NOX::Belos::Group::applyJacobianTranspose(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  return grpPtr->applyJacobianTranspose(input, result);
}

NOX::Abstract::Group::ReturnType
NOX::Belos::Group::applyJacobianInverse(NOX::Parameter::List& params,
					const NOX::Abstract::Vector& input,
					NOX::Abstract::Vector& result) const 
{
  // Create multivectors out of input, result
  NOX::Abstract::MultiVector *inputs = 
    input.createMultiVector(NULL, 0, NOX::DeepCopy);
  NOX::Abstract::MultiVector *results = 
    result.createMultiVector(NULL, 0, NOX::DeepCopy);
  
  // Call multivector version
  NOX::Abstract::Group::ReturnType res =
    applyJacobianInverseMultiVector(params, *inputs, *results);

  // Copy result
  result = (*results)[0];

  // Delete temporaries
  delete inputs;
  delete results;

  return res;
}

NOX::Abstract::Group::ReturnType 
NOX::Belos::Group::applyRightPreconditioning(
					 bool useTranspose,
					 NOX::Parameter::List& params,
					 const NOX::Abstract::Vector& input, 
					 NOX::Abstract::Vector& result) const
{
  return grpPtr->applyRightPreconditioning(useTranspose, params, input, 
					   result);
}

NOX::Abstract::Group::ReturnType
NOX::Belos::Group::applyJacobianMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  return grpPtr->applyJacobianMultiVector(input, result);
}

NOX::Abstract::Group::ReturnType
NOX::Belos::Group::applyJacobianTransposeMultiVector(
				     const NOX::Abstract::MultiVector& input,
				     NOX::Abstract::MultiVector& result) const 
{
  return grpPtr->applyJacobianTransposeMultiVector(input, result);
}

NOX::Abstract::Group::ReturnType
NOX::Belos::Group::applyJacobianInverseMultiVector(
				    NOX::Parameter::List& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const
{
  // Cast away const on input
  NOX::Abstract::MultiVector& nonConstInput = 
    const_cast<NOX::Abstract::MultiVector&>(input);

  // Create Belos Jacobian Operator
  NOX::Belos::JacobianOperator belosJacOp(*grpPtr);

  // Create Belos Preconditioner Operator
  NOX::Belos::PreconditionOperator belosPrecOp(*grpPtr, params);

  // Create Belos MultiVec's for input, result
  NOX::Belos::MultiVector belosInput(nonConstInput);
  NOX::Belos::MultiVector belosResult(result);

  // Create LinearProblemManager
  ::Belos::LinearProblemManager<double> belosProblem(&belosJacOp, 
						     &belosResult, 
						     &belosInput);
  belosProblem.SetRightPrec(&belosPrecOp);

  // Parse parameter list
  int maxits = params.getParameter("Max Iterations", 400);
  double tol = params.getParameter("Tolerance", 1.0e-6);
  int length = params.getParameter("Size of Krylov Subspace", 300);
  int numrestarts = params.getParameter("Number of Restarts", 20);
  int maxblocksize = params.getParameter("Maximum block size", 10);
  string method = params.getParameter("Belos Solver", "GMRES");
  int verbLevel = params.getParameter("Verbosity Level", 1);

  // Create status tests
  ::Belos::StatusTestMaxIters<double> test1(maxits);
  ::Belos::StatusTestMaxRestarts<double> test2(numrestarts);
  ::Belos::StatusTestCombo<double> belosTest(
					 ::Belos::StatusTestCombo<double>::OR,
					 test1, test2 );
  ::Belos::StatusTestResNorm<double> test3( tol );
  test3.DefineScaleForm(::Belos::StatusTestResNorm<double>::NormOfPrecInitRes, 
			::Belos::TwoNorm );
  belosTest.AddStatusTest( test3 );

  // Set block size
  int blocksize = input.numVectors();
  if (blocksize > maxblocksize)
    blocksize = maxblocksize;
  belosProblem.SetBlockSize(blocksize);

  // Create output manager
  ::Belos::OutputManager<double> belosOutputManager(myPID, verbLevel);

  if (method == "GMRES") {
    ::Belos::BlockGmres<double> belosGMRES(belosProblem, belosTest, 
					   belosOutputManager, length);
    belosGMRES.Solve();
  }
  else if (method == "CG") {
    ::Belos::BlockCG<double> belosCG(belosProblem, belosTest, 
				     belosOutputManager);
    belosCG.Solve();
  }
  else {
    cout << "ERROR: NOX::Belos::Group::applyJacobianInverseMultiVector" << endl
	 << "\"Belos Solver\" parameter \"" << method 
	 <<  "\" is invalid!" << endl;
    throw "NOX Error";
  }

  ::Belos::StatusType status = belosTest.GetStatus();
  if (status == ::Belos::Unconverged)
    return NOX::Abstract::Group::NotConverged;
  else if (status == ::Belos::Converged)
    return NOX::Abstract::Group::Ok;
  else 
    return NOX::Abstract::Group::Failed;
}

NOX::Abstract::Group::ReturnType
NOX::Belos::Group::applyRightPreconditioningMultiVector(
				   bool useTranspose,
				   NOX::Parameter::List& params,
				   const NOX::Abstract::MultiVector& input, 
				   NOX::Abstract::MultiVector& result) const
{
  return grpPtr->applyRightPreconditioningMultiVector(useTranspose, params, 
						      input, result);
}

bool
NOX::Belos::Group::isF() const 
{
  return grpPtr->isF();
}

bool
NOX::Belos::Group::isJacobian() const 
{
  return grpPtr->isJacobian();
}

bool
NOX::Belos::Group::isGradient() const 
{
  return grpPtr->isGradient();
}

bool
NOX::Belos::Group::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
NOX::Belos::Group::getX() const 
{
  return grpPtr->getX();
}

const NOX::Abstract::Vector&
NOX::Belos::Group::getF() const 
{
  return grpPtr->getF();
}

double
NOX::Belos::Group::getNormF() const 
{
  return grpPtr->getNormF();
}

const NOX::Abstract::Vector&
NOX::Belos::Group::getGradient() const 
{
  return grpPtr->getGradient();
}

const NOX::Abstract::Vector&
NOX::Belos::Group::getNewton() const 
{
  return *newtonVecPtr;
}

double
NOX::Belos::Group::getNormNewtonSolveResidual() const 
{
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Vector *residual = getF().clone(NOX::DeepCopy);
  
  status = applyJacobian(*newtonVecPtr, *residual);
  if (status != NOX::Abstract::Group::Ok) {
    cerr << "Error:  NOX::Belos::Group::getNormNewtonSolveResidual() -- applyJacobian failed!" << endl;
    throw "NOX Error";
  }

  residual->update(1.0, getF(), 1.0);
  double resid_norm = residual->norm();

  delete residual;

  return resid_norm;
}

void
NOX::Belos::Group::resetIsValid()
{
  isValidNewton = false;
}
