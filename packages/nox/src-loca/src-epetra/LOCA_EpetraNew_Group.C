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

#include "LOCA_EpetraNew_Group.H"	          // class definition

#include "LOCA_EpetraNew_Interface_Required.H"        // class data members
#include "NOX_Parameter_List.H"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Epetra_BorderedOp.H"
#include "LOCA_Epetra_HouseholderJacOp.H"

LOCA::EpetraNew::Group::Group(NOX::Parameter::List& printingParams, 
			      LOCA::EpetraNew::Interface::Required& i, 
			      NOX::Epetra::Vector& initialGuess,
			      const LOCA::ParameterVector& p) :
  NOX::EpetraNew::Group(printingParams, i, initialGuess),
  LOCA::Abstract::Group(),
  params(p),
  userInterfaceReq(i),
  userInterfaceTime(NULL),
  tmpVectorPtr2(0),
  scaleVecPtr(NULL)
{
}

LOCA::EpetraNew::Group::Group(NOX::Parameter::List& printingParams, 
			      LOCA::EpetraNew::Interface::Required& i, 
			      NOX::Epetra::Vector& initialGuess, 
			      NOX::EpetraNew::LinearSystem& linSys,
			      const LOCA::ParameterVector& p) :
  NOX::EpetraNew::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(),
  params(p),
  userInterfaceReq(i),
  userInterfaceTime(NULL),
  tmpVectorPtr2(0),
  scaleVecPtr(NULL)
{
}

LOCA::EpetraNew::Group::Group(NOX::Parameter::List& printingParams, 
			      LOCA::EpetraNew::Interface::TimeDependent& i, 
			      NOX::Epetra::Vector& initialGuess, 
			      NOX::EpetraNew::LinearSystem& linSys,
			      const LOCA::ParameterVector& p) :
  NOX::EpetraNew::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(),
  params(p),
  userInterfaceReq(i), 
  userInterfaceTime(&i),
  tmpVectorPtr2(0),
  scaleVecPtr(NULL)
{
}

LOCA::EpetraNew::Group::Group(const LOCA::EpetraNew::Group& source, 
			   NOX::CopyType type) :
  NOX::EpetraNew::Group(source, type),
  LOCA::Abstract::Group(source, type),
  params(source.params),
  userInterfaceReq(source.userInterfaceReq),
  userInterfaceTime(source.userInterfaceTime),
  tmpVectorPtr2(0),
  scaleVecPtr(NULL)
{
  if (source.scaleVecPtr != NULL)
    scaleVecPtr = source.scaleVecPtr->clone(NOX::DeepCopy);
}

LOCA::EpetraNew::Group::~Group() 
{
  delete tmpVectorPtr2;
  if (scaleVecPtr != NULL)
    delete scaleVecPtr;
}

NOX::Abstract::Group* 
LOCA::EpetraNew::Group::clone(NOX::CopyType type) const 
{
  return new LOCA::EpetraNew::Group(*this, type);
}

NOX::Abstract::Group& 
LOCA::EpetraNew::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::Abstract::Group& 
LOCA::EpetraNew::Group::operator=(const LOCA::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

NOX::Abstract::Group& 
LOCA::EpetraNew::Group::operator=(const NOX::EpetraNew::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::EpetraNew::Group& 
LOCA::EpetraNew::Group::operator=(const LOCA::EpetraNew::Group& source)
{
  params = source.params;
  NOX::EpetraNew::Group::operator=(source);
  LOCA::Abstract::Group::operator=(source);
  if (source.scaleVecPtr != NULL)
    scaleVecPtr = source.scaleVecPtr->clone(NOX::DeepCopy);
  return *this;
}

void 
LOCA::EpetraNew::Group::setParams(const LOCA::ParameterVector& p)
{
  resetIsValid();
  params = p;
}

void
LOCA::EpetraNew::Group::setParam(int paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::EpetraNew::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

void
LOCA::EpetraNew::Group::setParam(string paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::EpetraNew::Group::getParam(string paramID) const
{
  return params.getValue(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::EpetraNew::Group::computeF() 
{

  if (isF())
    return Abstract::Group::Ok;
  
  // Set the parameters prior to computing F
  userInterfaceReq.setParameters(params);
  
  return NOX::EpetraNew::Group::computeF();
}

NOX::Abstract::Group::ReturnType 
LOCA::EpetraNew::Group::computeJacobian() 
{

  if (isJacobian())
    return Abstract::Group::Ok;
  
  // Set the parameters prior to computing F
  userInterfaceReq.setParameters(params);

  return NOX::EpetraNew::Group::computeJacobian();
}

const LOCA::ParameterVector& 
LOCA::EpetraNew::Group::getParams() const 
{
  return params;
}

NOX::EpetraNew::Interface::Required& 
LOCA::EpetraNew::Group::getUserInterface()
{
  return userInterfaceReq;
}

void
LOCA::EpetraNew::Group::printSolution(const double conParam) const
{
  printSolution(xVector, conParam);
}

void
LOCA::EpetraNew::Group::printSolution(const NOX::Epetra::Vector& x_,
				      const double conParam) const
{
  userInterfaceReq.printSolution(x_.getEpetraVector(), conParam);
}

void
LOCA::EpetraNew::Group::printSolution(const NOX::Abstract::Vector& x_,
				      const double conParam) const
{
  printSolution(dynamic_cast<const NOX::Epetra::Vector&>(x_), conParam);
}

double
LOCA::EpetraNew::Group::computeScaledDotProduct(
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
   
void
LOCA::EpetraNew::Group::setScaleVector(const NOX::Abstract::Vector& s)
{
  if (scaleVecPtr != NULL)
    delete scaleVecPtr;

  scaleVecPtr = s.clone(NOX::DeepCopy);

  return;
}

NOX::Abstract::Group::ReturnType
LOCA::EpetraNew::Group::computeMassMatrix()
{
  if(userInterfaceTime != NULL)
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::BadDependency;
}

NOX::Abstract::Group::ReturnType
LOCA::EpetraNew::Group::applyMassMatrix(const NOX::Abstract::Vector& input,
                                        NOX::Abstract::Vector& result) const
{


  if(userInterfaceTime != NULL){
     const NOX::Epetra::Vector& epetraInput = 
       dynamic_cast<const NOX::Epetra::Vector&>(input);
     NOX::Epetra::Vector& epetraResult =
       dynamic_cast<NOX::Epetra::Vector&>(result);
     dynamic_cast<LOCA::EpetraNew::Interface::TimeDependent&>(userInterfaceReq).applyMassMatrix(epetraInput,epetraResult);
  }
  else
    return NOX::Abstract::Group::BadDependency;
}

bool 
LOCA::EpetraNew::Group::isMassMatrix()
{
  return true;
}

NOX::Abstract::Group::ReturnType 
LOCA::EpetraNew::Group::applyHouseholderJacobianInverse(
					  NOX::Parameter::List& p,
					  const NOX::Abstract::Vector& cf,
					  const NOX::Abstract::Vector& cdfdp,
					  const NOX::Abstract::Vector& cux,
					  double up, double beta,
					  NOX::Abstract::Vector& result_x,
					  double& result_p) const
{
  // Get non-const input
  NOX::Abstract::Vector& f = const_cast<NOX::Abstract::Vector&>(cf);
  NOX::Abstract::Vector& dfdp = const_cast<NOX::Abstract::Vector&>(cdfdp);
  NOX::Abstract::Vector& ux = const_cast<NOX::Abstract::Vector&>(cux);

  // cast vectors to nox epetra vectors
  NOX::Epetra::Vector& nox_epetra_f = 
    dynamic_cast<NOX::Epetra::Vector&>(f);
  NOX::Epetra::Vector& nox_epetra_dfdp = 
    dynamic_cast<NOX::Epetra::Vector&>(dfdp);
  NOX::Epetra::Vector& nox_epetra_ux = 
    dynamic_cast<NOX::Epetra::Vector&>(ux);
  NOX::Epetra::Vector& nox_epetra_result_x = 
    dynamic_cast<NOX::Epetra::Vector&>(result_x);
  
  // Get underlying epetra vectors
  Epetra_Vector& epetra_f = nox_epetra_f.getEpetraVector();
  Epetra_Vector& epetra_dfdp = nox_epetra_dfdp.getEpetraVector();
  Epetra_Vector& epetra_ux = nox_epetra_ux.getEpetraVector();
  Epetra_Vector& epetra_result_x = nox_epetra_result_x.getEpetraVector();

  // Get Jacobian operator
  const Epetra_Operator& jac = 
    sharedLinearSystem.getObject(this).getJacobianOperator();

  // Build Householder projection-Jacobian operator
  LOCA::Epetra::HouseholderJacOp houseJac(jac, epetra_dfdp, epetra_ux,
					  up, beta);

  // Initialize operator
  houseJac.init(epetra_result_x);

  // Set Jacobian operator for solve
  sharedLinearSystem.getObject(this).setJacobianOperatorForSolve(houseJac);

  // compute preconditioner if necessary
  bool reusePrec = 
    sharedLinearSystem.getObject(this).checkPreconditionerReuse();

  if (!isPreconditioner()  && !reusePrec ) {
    sharedLinearSystem.getObject(this).destroyPreconditioner();
    sharedLinearSystem.getObject(this)
      .createPreconditioner(xVector.getEpetraVector(), p, false);
    isValidPreconditioner = true;
  }

  bool status = 
    sharedLinearSystem.getObject(this).applyJacobianInverse(
							 p, 
							 nox_epetra_f, 
							 nox_epetra_result_x);

  // Apply Householder transformation to result
  result_p = 0.0;
  houseJac.applyHouse(epetra_result_x, result_p);

  // Finalize operator
  houseJac.finish();

  if (status) 
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::NotConverged;
}

void
LOCA::EpetraNew::Group::scaleVector(NOX::Abstract::Vector& x) const
{
  if (scaleVecPtr == NULL)
    x.scale(1.0 / sqrt(static_cast<double>(x.length())));
  else 
    x.scale(*scaleVecPtr);
}

void
LOCA::EpetraNew::Group::projectToDraw(const NOX::Abstract::Vector& x,
				   double *px) const
{
  const NOX::Epetra::Vector& ex = 
    dynamic_cast<const NOX::Epetra::Vector&>(x);
  userInterfaceReq.projectToDraw(ex, px);
}

int
LOCA::EpetraNew::Group::projectToDrawDimension() const
{
  return userInterfaceReq.projectToDrawDimension();
}

NOX::Abstract::Group::ReturnType 
LOCA::EpetraNew::Group::augmentJacobianForHomotopy(double conParamValue)
{

  //Allocate temporary vectors if not aready done
  if (tmpVectorPtr == 0)
    tmpVectorPtr = new Epetra_Vector(xVector.getEpetraVector());
  if (tmpVectorPtr2 == 0)
    tmpVectorPtr2 = new Epetra_Vector(xVector.getEpetraVector());

  tmpVectorPtr2->PutScalar(1.0-conParamValue);

  // See if it is an Epetra_CrsMatrix
  const Epetra_CrsMatrix* constTestCrs = 0;
  Epetra_CrsMatrix* testCrs = 0;
  constTestCrs = dynamic_cast<const Epetra_CrsMatrix*>
    (&(sharedLinearSystem.getObject(this).getJacobianOperator()));
  if (constTestCrs != 0) {
    testCrs = const_cast<Epetra_CrsMatrix*>(constTestCrs);
    testCrs->Scale(conParamValue);
    testCrs->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testCrs->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;

  }

  // See if it is an Epetra_VbrMatrix
  const Epetra_VbrMatrix* constTestVbr = 0;
  Epetra_VbrMatrix* testVbr = 0;
  constTestVbr = dynamic_cast<const Epetra_VbrMatrix*>
    (&(sharedLinearSystem.getObject(this).getJacobianOperator()));
  if (testVbr != 0) {
    testVbr = const_cast<Epetra_VbrMatrix*>(constTestVbr);
    testVbr->Scale(conParamValue);
    testVbr->ExtractDiagonalCopy(*tmpVectorPtr);
    tmpVectorPtr->Update(1.0, *tmpVectorPtr2, 1.0);
    testVbr->ReplaceDiagonalValues(*tmpVectorPtr);
    return LOCA::Abstract::Group::Ok;
  }

  // Otherwise this alg won't work -- return NotDefined

  return LOCA::Abstract::Group::NotDefined;
}

// NOX::Abstract::Group::ReturnType 
// LOCA::EpetraNew::Group::applyBorderedJacobianInverse(
// 				     bool trans,
// 				     NOX::Parameter::List& p,
// 				     const NOX::Abstract::Vector& ca,
// 				     const NOX::Abstract::Vector& cb,
// 				     const NOX::Abstract::Vector& cvInput,
// 				     double sInput,
// 				     NOX::Abstract::Vector& vResult,
// 				     double& sResult) const
// {
//   // Get non-const input
//   NOX::Abstract::Vector& vInput = const_cast<NOX::Abstract::Vector&>(cvInput);
//   NOX::Abstract::Vector& a = const_cast<NOX::Abstract::Vector&>(ca);
//   NOX::Abstract::Vector& b = const_cast<NOX::Abstract::Vector&>(cb);

//   // cast vectors to nox epetra vectors
//   NOX::Epetra::Vector& nox_epetra_vInput = 
//     dynamic_cast<NOX::Epetra::Vector&>(vInput);
//   NOX::Epetra::Vector& nox_epetra_a = 
//     dynamic_cast<NOX::Epetra::Vector&>(a);
//   NOX::Epetra::Vector& nox_epetra_b = 
//     dynamic_cast<NOX::Epetra::Vector&>(b);
//   NOX::Epetra::Vector& nox_epetra_vResult = 
//     dynamic_cast<NOX::Epetra::Vector&>(vResult);
  
//   // Get underlying epetra vectors
//   Epetra_Vector& epetra_vInput = nox_epetra_vInput.getEpetraVector();
//   Epetra_Vector& epetra_a = nox_epetra_a.getEpetraVector();
//   Epetra_Vector& epetra_b = nox_epetra_b.getEpetraVector();
//   Epetra_Vector& epetra_vResult = nox_epetra_vResult.getEpetraVector();

//   // Get Jacobian, preconditioner operators
//   const Epetra_Operator& cjac = 
//     sharedLinearSystem.getObject(this).getJacobianOperator();
//   const Epetra_Operator& cprec = 
//     sharedLinearSystem.getObject(this).getGeneratedPrecOperator();
//   Epetra_Operator& jac = const_cast<Epetra_Operator&>(cjac);
//   Epetra_Operator& prec = const_cast<Epetra_Operator&>(cprec);

//   // Build bordered matrix-free Jacobian, preconditioning operator
//   LOCA::Epetra::BorderedOp extended_jac(jac, epetra_a, epetra_b);
//   LOCA::Epetra::BorderedOp extended_prec(prec, epetra_a, epetra_b);
//   extended_jac.SetUseTranspose(trans);
//   extended_prec.SetUseTranspose(trans);

//   // Build extended epetra vectors
//   Epetra_Vector *epetra_extended_input = 
//     extended_jac.buildEpetraExtendedVec(epetra_vInput, sInput, true);
//   Epetra_Vector *epetra_extended_result = 
//     extended_jac.buildEpetraExtendedVec(epetra_vResult, 0.0, false);

//   // Build extended NOX::Epetra vectors
//   NOX::Epetra::Vector nox_epetra_extended_input(*epetra_extended_input,
// 						NOX::DeepCopy,
// 						true);
//   NOX::Epetra::Vector nox_epetra_extended_result(*epetra_extended_result,
// 						 NOX::DeepCopy,
// 						 true);

//   bool status = 
//     sharedLinearSystem.getObject(this).applyJacobianInverse(
// 						   p, 
// 						   nox_epetra_extended_input, 
// 						   nox_epetra_extended_result,
// 						   &extended_jac,
// 						   &extended_prec);

//   extended_jac.setEpetraExtendedVec(epetra_vResult, sResult, 
// 				    *epetra_extended_result);
  
//   delete epetra_extended_input;
//   delete epetra_extended_result;

//   if (status) 
//     return NOX::Abstract::Group::Ok;
//   else
//     return NOX::Abstract::Group::NotConverged;
// }
