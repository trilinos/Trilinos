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

#include "LOCA_Epetra_Group.H"	          // class definition

#include "LOCA_Epetra_Interface.H"        // class data members
#include "NOX_Parameter_List.H"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "NOX_Epetra_Scaling.H"
#include "NOX_Epetra_SharedOperator.H"
#include "Epetra_Vector.h"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Epetra_BorderedOp.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h" 
#include "Epetra_Vector.h"

// External include files - linking to Aztec00 and Epetra in Trilinos
#include "AztecOO.h"
#include "AztecOO_Operator.h"
#include "AztecOO_StatusTest.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_StatusTestMaxIters.h"
#include "AztecOO_StatusTestResNorm.h"

#include "Ifpack_CrsRiluk.h"

LOCA::Epetra::Group::Group(NOX::Parameter::List& printParams,
			   NOX::Parameter::List& par, 
			   LOCA::Epetra::Interface& i, 
			   const LOCA::ParameterVector& p, 
			   NOX::Epetra::Vector& x, 
			   Epetra_Operator& J) :
  NOX::Epetra::Group(printParams, par, i, x, J),
  LOCA::Abstract::Group(par),
  params(p),
  userInterface(i),
  tmpVectorPtr2(0),
  scaleVecPtr(NULL)
{
}

LOCA::Epetra::Group::Group(NOX::Parameter::List& printParams,
			   NOX::Parameter::List& par, 
			   LOCA::Epetra::Interface& i, 
			   const LOCA::ParameterVector& p, 
			   NOX::Epetra::Vector& x, 
			   Epetra_Operator& J, 
			   Epetra_Operator& M) :
  NOX::Epetra::Group(printParams, par, i, x, J, M),
  LOCA::Abstract::Group(par),
  params(p),
  userInterface(i),
  tmpVectorPtr2(0),
  scaleVecPtr(NULL)
{
}

LOCA::Epetra::Group::Group(const LOCA::Epetra::Group& source, 
			   NOX::CopyType type) :
  NOX::Epetra::Group(source, type),
  LOCA::Abstract::Group(source, type),
  params(source.params),
  userInterface(source.userInterface),
  tmpVectorPtr2(0),
  scaleVecPtr(NULL)
{
  if (source.scaleVecPtr != NULL)
    scaleVecPtr = source.scaleVecPtr->clone(NOX::DeepCopy);
}

LOCA::Epetra::Group::~Group() 
{
  delete tmpVectorPtr2;
  if (scaleVecPtr != NULL)
    delete scaleVecPtr;
}

NOX::Abstract::Group* 
LOCA::Epetra::Group::clone(NOX::CopyType type) const 
{
  return new Group(*this, type);
}

NOX::Abstract::Group& 
LOCA::Epetra::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::Abstract::Group& 
LOCA::Epetra::Group::operator=(const LOCA::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::Epetra::Group& 
LOCA::Epetra::Group::operator=(const LOCA::Epetra::Group& source)
{
  params = source.params;
  NOX::Epetra::Group::operator=(source);
  LOCA::Abstract::Group::operator=(source);
  if (source.scaleVecPtr != NULL)
    scaleVecPtr = source.scaleVecPtr->clone(NOX::DeepCopy);
  return *this;
}

void 
LOCA::Epetra::Group::setParams(const LOCA::ParameterVector& p)
{
  resetIsValid();
  params = p;
}

void
LOCA::Epetra::Group::setParam(int paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::Epetra::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

void
LOCA::Epetra::Group::setParam(string paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::Epetra::Group::getParam(string paramID) const
{
  return params.getValue(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeF() 
{

  if (isF())
    return Abstract::Group::Ok;
  
  // Set the parameters prior to computing F
  userInterface.setParameters(params);
  
  return NOX::Epetra::Group::computeF();
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::computeJacobian() 
{

  if (isJacobian())
    return Abstract::Group::Ok;
  
  // Set the parameters prior to computing F
  userInterface.setParameters(params);

  return NOX::Epetra::Group::computeJacobian();
}

const LOCA::ParameterVector& 
LOCA::Epetra::Group::getParams() const 
{
  return params;
}

NOX::Epetra::Interface& 
LOCA::Epetra::Group::getUserInterface()
{
  return userInterface;
}

void
LOCA::Epetra::Group::printSolution(const double conParam) const
{
  printSolution(xVector, conParam);
}

void
LOCA::Epetra::Group::printSolution(const NOX::Epetra::Vector& x_,
                                   const double conParam) const
{
  userInterface.printSolution(x_.getEpetraVector(), conParam);
}

void
LOCA::Epetra::Group::printSolution(const NOX::Abstract::Vector& x_,
                                   const double conParam) const
{
  printSolution(dynamic_cast<const NOX::Epetra::Vector&>(x_), conParam);
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::augmentJacobianForHomotopy(double conParamValue)
{

  //Allocate temporary vectors if not aready done
  if (tmpVectorPtr == 0)
    tmpVectorPtr = new Epetra_Vector(xVector.getEpetraVector());
  if (tmpVectorPtr2 == 0)
    tmpVectorPtr2 = new Epetra_Vector(xVector.getEpetraVector());

  tmpVectorPtr2->PutScalar(1.0-conParamValue);

  // See if it is an Epetra_CrsMatrix
  Epetra_CrsMatrix* testCrs = 0;
  testCrs = dynamic_cast<Epetra_CrsMatrix*>
            (&(sharedJacobian.getOperator(this)));
  if (testCrs != 0) {

    testCrs->Scale(conParamValue);
    testCrs->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testCrs->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;

  }

  // See if it is an Epetra_VbrMatrix
  Epetra_VbrMatrix* testVbr = 0;
  testVbr = dynamic_cast<Epetra_VbrMatrix*>(&(sharedJacobian.getOperator(this)));
  if (testVbr != 0) {
    
    testVbr->Scale(conParamValue);
    testVbr->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testVbr->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;
  }

  // Otherwise this alg won't work
  LOCA::ErrorCheck::throwError(
    "LOCA::Epetra::Group::augmentJacobianForHomotopy()",
    "the Jacobian must be either an Epetra_CrsMatrix or an Epetra_VbrMatrix!");

  return LOCA::Abstract::Group::Ok;
}

double
LOCA::Epetra::Group::computeScaledDotProduct(
				       const NOX::Abstract::Vector& a,
				       const NOX::Abstract::Vector& b) const
{
  if (scaleVecPtr == NULL)
    return a.dot(b) / a.length();
  else {
    NOX::Abstract::Vector* as = a.clone(NOX::DeepCopy);
    NOX::Abstract::Vector* bs = b.clone(NOX::DeepCopy);
    double d;

    as->scale(*scaleVecPtr);
    bs->scale(*scaleVecPtr);
    d = as->dot(*bs);

    delete as;
    delete bs;

    return d;
  }
}
    
NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::applyBorderedJacobianInverse(bool trans,
				     NOX::Parameter::List& p,
				     const NOX::Abstract::Vector& ca,
				     const NOX::Abstract::Vector& cb,
				     const NOX::Abstract::Vector& cvInput,
				     double sInput,
				     NOX::Abstract::Vector& vResult,
				     double& sResult) const
{
  // Get non-const input
  NOX::Abstract::Vector& vInput = const_cast<NOX::Abstract::Vector&>(cvInput);
  NOX::Abstract::Vector& a = const_cast<NOX::Abstract::Vector&>(ca);
  NOX::Abstract::Vector& b = const_cast<NOX::Abstract::Vector&>(cb);

  // cast vectors to nox epetra vectors
  NOX::Epetra::Vector& nox_epetra_vInput = 
    dynamic_cast<NOX::Epetra::Vector&>(vInput);
  NOX::Epetra::Vector& nox_epetra_a = 
    dynamic_cast<NOX::Epetra::Vector&>(a);
  NOX::Epetra::Vector& nox_epetra_b = 
    dynamic_cast<NOX::Epetra::Vector&>(b);
  NOX::Epetra::Vector& nox_epetra_vResult = 
    dynamic_cast<NOX::Epetra::Vector&>(vResult);
  
  // Get underlying epetra vectors
  Epetra_Vector& epetra_vInput = nox_epetra_vInput.getEpetraVector();
  Epetra_Vector& epetra_a = nox_epetra_a.getEpetraVector();
  Epetra_Vector& epetra_b = nox_epetra_b.getEpetraVector();
  Epetra_Vector& epetra_vResult = nox_epetra_vResult.getEpetraVector();

  // Create preconditioner
  createIfpackPreconditioner(p);

  // Get Jacobian, preconditioner operators
  //const Epetra_Operator& cjac = sharedJacobianPtr->getOperator();
  const Epetra_Operator& cjac = sharedJacobian.getOperator();
  Epetra_Operator& jac = const_cast<Epetra_Operator&>(cjac);
  Epetra_Operator& prec = *ifpackPreconditioner;

  // Build bordered matrix-free Jacobian, preconditioning operator
  LOCA::Epetra::BorderedOp extended_jac(jac, epetra_a, epetra_b);
  LOCA::Epetra::BorderedOp extended_prec(prec, epetra_a, epetra_b);
  extended_jac.SetUseTranspose(trans);
  extended_prec.SetUseTranspose(trans);

  // Build extended epetra vectors
  Epetra_Vector *epetra_extended_input = 
    extended_jac.buildEpetraExtendedVec(epetra_vInput, sInput, true);
  Epetra_Vector *epetra_extended_result = 
    extended_jac.buildEpetraExtendedVec(epetra_vResult, 0.0, false);

  // Create Epetra linear problem object for the linear solve
  Epetra_LinearProblem Problem(&extended_jac, 
  			       epetra_extended_result, 
			       epetra_extended_input);

  // Set the default Problem parameters to "hard" (this sets Aztec defaults
  // during the AztecOO instantiation)
  Problem.SetPDL(hard);

  // Create the solver. 
  AztecOO borderedAztecSolver(Problem);

  // Set specific Aztec parameters requested by NOX
  setBorderedAztecOptions(p, borderedAztecSolver);

  // Set preconditioner
  borderedAztecSolver.SetPrecOperator(&extended_prec);

  // Get linear solver convergence parameters
  int maxit = p.getParameter("Max Iterations", 400);
  double tol = p.getParameter("Tolerance", 1.0e-6);
  
  // Solve linear problem
  int aztecStatus = borderedAztecSolver.Iterate(maxit, tol);

  // Set the output parameters in the "Output" sublist
  NOX::Parameter::List& outputList = p.sublist("Output");
  int prevLinIters = 
    outputList.getParameter("Total Number of Linear Iterations", 0);
  int curLinIters = borderedAztecSolver.NumIters();
  double achievedTol = borderedAztecSolver.ScaledResidual();

  outputList.setParameter("Number of Linear Iterations", curLinIters);
  outputList.setParameter("Total Number of Linear Iterations", 
			  (prevLinIters + curLinIters));
  outputList.setParameter("Achieved Tolerance", achievedTol);

  extended_jac.setEpetraExtendedVec(epetra_vResult, sResult, 
				    *epetra_extended_result);
  
  delete epetra_extended_input;
  delete epetra_extended_result;

  // Destroy preconditioner
  destroyPreconditioner();

  if (aztecStatus != 0) 
    return NOX::Abstract::Group::NotConverged;
  
  return NOX::Abstract::Group::Ok;
}

void
LOCA::Epetra::Group::setScaleVector(const NOX::Abstract::Vector& s)
{
  if (scaleVecPtr != NULL)
    delete scaleVecPtr;

  scaleVecPtr = s.clone(NOX::DeepCopy);

  return;
}

void 
LOCA::Epetra::Group::setBorderedAztecOptions(const NOX::Parameter::List& p, 
					     AztecOO& aztec) const
{
  // Set the Aztec Solver
  string linearSolver = p.getParameter("Aztec Solver", "GMRES");
  if (linearSolver == "CG")
    aztec.SetAztecOption(AZ_solver, AZ_cg);
  else if (linearSolver == "GMRES")
    aztec.SetAztecOption(AZ_solver, AZ_gmres);
  else if (linearSolver == "CGS")
    aztec.SetAztecOption(AZ_solver, AZ_cgs);
  else if (linearSolver == "TFQMR")
    aztec.SetAztecOption(AZ_solver, AZ_tfqmr);
  else if (linearSolver == "BiCGStab")
    aztec.SetAztecOption(AZ_solver, AZ_bicgstab);
  else if (linearSolver == "LU")
    aztec.SetAztecOption(AZ_solver, AZ_lu);
  else {
    cout << "ERROR: LOCA::Epetra::Group::setAztecOptions" 
	 << endl
	 << "\"Aztec Solver\" parameter \"" << linearSolver 
	 <<  "\" is invalid!" << endl;
    throw "NOX Error";
  }
    
  // Gram-Schmidt orthogonalization procedure
  string orthog = p.getParameter("Orthogonalization", "Classical");
  if (orthog == "Classical") 
    aztec.SetAztecOption(AZ_orthog, AZ_classic);
  else if (orthog == "Modified")
    aztec.SetAztecOption(AZ_orthog, AZ_modified);
  else {
    cout << "ERROR: LOCA::Epetra::Group::setAztecOptions" 
	 << endl
	 << "\"Orthogonalization\" parameter \"" << orthog
	 << "\" is invalid!" << endl;
    throw "NOX Error";
  }

  // Size of the krylov space
  aztec.SetAztecOption(AZ_kspace, 
		       p.getParameter("Size of Krylov Subspace", 300));

  // Convergence criteria to use in the linear solver
  string convCriteria = p.getParameter("Convergence Criteria", "r0");
  if (convCriteria == "r0") 
    aztec.SetAztecOption(AZ_conv, AZ_r0);
  else if (convCriteria == "rhs")
    aztec.SetAztecOption(AZ_conv, AZ_rhs);
  else if (convCriteria == "Anorm")
    aztec.SetAztecOption(AZ_conv, AZ_Anorm);
  else if (convCriteria == "no scaling")
    aztec.SetAztecOption(AZ_conv, AZ_noscaled);
  else if (convCriteria == "sol")
    aztec.SetAztecOption(AZ_conv, AZ_sol);
  else {
    cout << "ERROR: LOCA::Epetra::Group::setAztecOptions" 
	 << endl
	 << "\"Convergence Criteria\" parameter \"" << convCriteria
	 << "\" is invalid!" << endl;
    throw "NOX Error";
  }

  // Set the ill-conditioning threshold for the upper hessenberg matrix
  if (p.isParameter("Ill-Conditioning Threshold")) {
    aztec.SetAztecParam(AZ_ill_cond_thresh, 
			p.getParameter("Ill-Conditioning Threshold", 1.0e+11));
  }

  // Frequency of linear solve residual output
  aztec.SetAztecOption(AZ_output, 
		       p.getParameter("Output Frequency", AZ_last));

  return;
}

