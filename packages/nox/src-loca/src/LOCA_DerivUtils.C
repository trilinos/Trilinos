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


// Class Header 
#include "LOCA_DerivUtils.H"

#include "LOCA_Continuation_AbstractGroup.H" 
#include "LOCA_Bifurcation_HopfBord_AbstractGroup.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_MultiVector.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Common.H"  // For fabs function
#include "LOCA_ErrorCheck.H"

LOCA::DerivUtils::DerivUtils(double perturb) :
  perturb(perturb)
{
  // assert (perturb>0.0);
}

LOCA::DerivUtils::DerivUtils(const DerivUtils& source) :
  perturb(source.perturb)
{

}

LOCA::DerivUtils::~DerivUtils()
{
 
}

LOCA::DerivUtils* 
LOCA::DerivUtils::clone(NOX::CopyType type) const
{
  return new DerivUtils(*this);  //Call Copy Constructor
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDfDp(LOCA::Continuation::AbstractGroup& grp,
			      int param_id, 
			      NOX::Abstract::Vector& result) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDfDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Compute base residual F
  if (!grp.isF()) {
    finalStatus = grp.computeF();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  // Copy original residual vector
  NOX::Abstract::Vector *Fvec = grp.getF().clone(NOX::DeepCopy);
  
  // Perturb single parameter in this group, and return perturbation, eps
  double param;
  double eps = perturbParam(grp, param, param_id);

  // Compute perturbed residual
  status = grp.computeF(); 
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Difference perturbed and base vector 
  result.update(1.0, grp.getF(), -1.0, *Fvec, 0.0);
  result.scale(1.0/eps);

  delete Fvec;

  // Restore original parameter value
  grp.setParam(param_id, param);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDfDp(LOCA::Continuation::AbstractGroup& grp,
			      const vector<int>& param_ids,
			      NOX::Abstract::MultiVector& result,
			      bool isValidF) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDfDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Views of f, df/dp
  NOX::Abstract::Vector *f = &result[0];
  NOX::Abstract::Vector *dfdp = NULL;

  // Compute base residual F
  if (!isValidF) {
    finalStatus = grp.computeF();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
    *f = grp.getF();
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;
  
  double param;
  double eps;

  // Loop over each parameter
  for (unsigned int i=0; i<param_ids.size(); i++) {

    // Perturb single parameter in this group, and return perturbation, eps
    eps = perturbParam(grp, param, param_ids[i]);

    // Compute perturbed residual
    status = grp.computeF(); 
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

    // Difference perturbed and base vector 
    dfdp = &result[i+1];
    dfdp->update(1.0, grp.getF(), -1.0, *f, 0.0);
    dfdp->scale(1.0/eps);

    // Restore original parameter value
    grp.setParam(param_ids[i], param);

  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDJnDp(LOCA::Continuation::AbstractGroup& grp,
			       const NOX::Abstract::Vector& nullVector, 
			       const int param_id,
			       NOX::Abstract::Vector& result) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDJnDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Allocate base Jn vector and fill with J time n
  NOX::Abstract::Vector *baseJnVectorPtr = nullVector.clone(NOX::ShapeCopy);

  if (!grp.isJacobian()) {
    finalStatus = grp.computeJacobian();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  status = grp.applyJacobian(nullVector, *baseJnVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Now that Jn is known, call other routine
  status = computeDJnDp(grp, nullVector, param_id, *baseJnVectorPtr, result);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  delete baseJnVectorPtr; 

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDJnDp(LOCA::Continuation::AbstractGroup& grp,
			       const NOX::Abstract::Vector& nullVector, 
			       const int param_id,
			       const NOX::Abstract::Vector& JnVector,
			       NOX::Abstract::Vector& result) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDJnDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Perturb single parameter in this group, and return perturbation
  double param;
  double eps = perturbParam(grp, param, param_id);

  // Fill perturbed Jn vector
  finalStatus = grp.computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  result.init(0.0);
  status = grp.applyJacobian(nullVector, result);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Difference perturbed and base vector 
  result.update(-1.0, JnVector, 1.0);
  result.scale(1.0/eps);

  // Restore original parameter value
  grp.setParam(param_id, param);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDJnDxa(LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector& nullVector,
				const NOX::Abstract::Vector& aVector,
				NOX::Abstract::Vector& result) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDJnDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Allocate base Jn vector and fill with J times n
  NOX::Abstract::Vector *baseJnVectorPtr = nullVector.clone(NOX::ShapeCopy);
  
  if (!grp.isJacobian()) {
    finalStatus = grp.computeJacobian();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  status = grp.applyJacobian(nullVector, *baseJnVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Now that Jn is known, call other routine
  status = computeDJnDxa(grp, nullVector, aVector, *baseJnVectorPtr, result);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  delete baseJnVectorPtr; 

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDJnDxa(LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector& nullVector,
				const NOX::Abstract::Vector& aVector,
				const NOX::Abstract::Vector& JnVector,
				NOX::Abstract::Vector& result) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDJnDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Copy original solution vector
  NOX::Abstract::Vector *Xvec = grp.getX().clone(NOX::DeepCopy);

  // Perturb solution vector in direction of aVector, return perturbation
  double eps = perturbXVec(grp, *Xvec, aVector);

  // Fill perturbed Jn vector
  finalStatus = grp.computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  
  result.init(0.0);
  status = grp.applyJacobian(nullVector, result);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Difference perturbed and base vector 
  result.update(-1.0, JnVector, 1.0);
  result.scale(1.0/eps);

  // Restore original solution vector
  grp.setX(*Xvec);

  delete Xvec;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDJnDxa(LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector& nullVector,
				const NOX::Abstract::MultiVector& aVector,
				NOX::Abstract::MultiVector& result) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDJnDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Allocate base Jn vector and fill with J times n
  NOX::Abstract::Vector *baseJnVectorPtr = nullVector.clone(NOX::ShapeCopy);
  
  if (!grp.isJacobian()) {
    finalStatus = grp.computeJacobian();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  status = grp.applyJacobian(nullVector, *baseJnVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Now that Jn is known, call other routine
  status = computeDJnDxa(grp, nullVector, aVector, *baseJnVectorPtr, result);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  delete baseJnVectorPtr; 

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDJnDxa(LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector& nullVector,
				const NOX::Abstract::MultiVector& aVector,
				const NOX::Abstract::Vector& JnVector,
				NOX::Abstract::MultiVector& result) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDJnDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Copy original solution vector
  NOX::Abstract::Vector *Xvec = grp.getX().clone(NOX::DeepCopy);

  // Loop over each column of multivector
  for (int i=0; i<aVector.numVectors(); i++) {

    // Perturb solution vector in direction of aVector, return perturbation
    double eps = perturbXVec(grp, *Xvec, aVector[i]);

    // Fill perturbed Jn vector
    finalStatus = grp.computeJacobian();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
    
    status = grp.applyJacobian(nullVector, result[i]);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

    // Difference perturbed and base vector 
    result[i].update(-1.0, JnVector, 1.0);
    result[i].scale(1.0/eps);

  }
  
  // Restore original solution vector
  grp.setX(*Xvec);

  delete Xvec;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDCeDp(
			   LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
			   const NOX::Abstract::Vector& yVector,
			   const NOX::Abstract::Vector& zVector,
			   double w,
			   const int param_id, 
			   NOX::Abstract::Vector& result_real,
			   NOX::Abstract::Vector& result_imag) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDCeDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Allocate base Ce vectors and fill with C times e
  NOX::Abstract::Vector *baseCeRealVectorPtr = yVector.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *baseCeImagVectorPtr = yVector.clone(NOX::ShapeCopy);

  if (!grp.isJacobian()) {
    finalStatus = grp.computeJacobian();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  if (!grp.isMassMatrix()) {
    status = grp.computeMassMatrix();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  else
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(NOX::Abstract::Group::Ok, 
						   finalStatus,
						   callingFunction);

  status = grp.applyComplex(yVector, zVector, w, *baseCeRealVectorPtr, 
			    *baseCeImagVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Now that Ce is known, call other routine
  status = 
    computeDCeDp(grp, yVector, zVector, w, param_id, *baseCeRealVectorPtr, 
		 *baseCeImagVectorPtr, result_real, result_imag);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  delete baseCeRealVectorPtr;
  delete baseCeImagVectorPtr;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDCeDp(
			   LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
			   const NOX::Abstract::Vector& yVector,
			   const NOX::Abstract::Vector& zVector,
			   double w,
			   const int param_id, 
			   const NOX::Abstract::Vector& Ce_real,
			   const NOX::Abstract::Vector& Ce_imag,
			   NOX::Abstract::Vector& result_real,
			   NOX::Abstract::Vector& result_imag) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDCeDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;


  // Perturb single parameter in this group, and return perturbation
  double param;
  double eps = perturbParam(grp, param, param_id);

  // Compute perturbed Ce vectors
  finalStatus = grp.computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  status = grp.computeMassMatrix();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  status = 
    grp.applyComplex(yVector, zVector, w, result_real,
				  result_imag);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Difference perturbed and base vector and return approximate derivative
  result_real.update(-1.0, Ce_real, 1.0); result_real.scale(1.0/eps);
  result_imag.update(-1.0, Ce_imag, 1.0); result_real.scale(1.0/eps);

  // Restore original parameter value
  grp.setParam(param_id, param);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDCeDxa(
			    LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
			    const NOX::Abstract::Vector& yVector,
			    const NOX::Abstract::Vector& zVector,
			    double w,
			    const NOX::Abstract::Vector& aVector,
			    NOX::Abstract::Vector& result_real,
			    NOX::Abstract::Vector& result_imag) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDCeDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Allocate base Ce vectors and fill with C times e
  NOX::Abstract::Vector *baseCeRealVectorPtr = yVector.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *baseCeImagVectorPtr = yVector.clone(NOX::ShapeCopy);

  if (!grp.isJacobian()) {
    finalStatus = grp.computeJacobian();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  if (!grp.isMassMatrix()) {
    status = grp.computeMassMatrix();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  else
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(NOX::Abstract::Group::Ok, 
						   finalStatus,
						   callingFunction);
  
  status = grp.applyComplex(yVector, zVector, w, *baseCeRealVectorPtr, 
			    *baseCeImagVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Now that Ce is known, call other routine
  status = 
    computeDCeDxa(grp, yVector, zVector, w, aVector, *baseCeRealVectorPtr, 
		  *baseCeImagVectorPtr, result_real, result_imag);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  delete baseCeRealVectorPtr; 
  delete baseCeImagVectorPtr; 

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDCeDxa(
			    LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
			    const NOX::Abstract::Vector& yVector,
			    const NOX::Abstract::Vector& zVector,
			    double w,
			    const NOX::Abstract::Vector& aVector,
			    const NOX::Abstract::Vector& Ce_real,
			    const NOX::Abstract::Vector& Ce_imag,
			    NOX::Abstract::Vector& result_real,
			    NOX::Abstract::Vector& result_imag) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDCeDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Copy original solution vector
  NOX::Abstract::Vector *Xvec = grp.getX().clone(NOX::DeepCopy);

  // Perturb solution vector in direction of aVector, return perturbation
  double eps = perturbXVec(grp, *Xvec, aVector);

  // Compute perturbed Ce vectors
  finalStatus = grp.computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  status = grp.computeMassMatrix();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  status = 
    grp.applyComplex(yVector, zVector, w, result_real,
				  result_imag);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Difference perturbed and base vector and return approximate derivative
  result_real.update(-1.0, Ce_real, 1.0); result_real.scale(1.0/eps);
  result_imag.update(-1.0, Ce_imag, 1.0); result_imag.scale(1.0/eps);

  // Restore original solution vector
  grp.setX(*Xvec);

  delete Xvec;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDCeDxa(
			    LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
			    const NOX::Abstract::Vector& yVector,
			    const NOX::Abstract::Vector& zVector,
			    double w,
			    const NOX::Abstract::MultiVector& aVector,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDCeDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Allocate base Ce vectors and fill with C times e
  NOX::Abstract::Vector *baseCeRealVectorPtr = yVector.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *baseCeImagVectorPtr = yVector.clone(NOX::ShapeCopy);

  if (!grp.isJacobian()) {
    finalStatus = grp.computeJacobian();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  if (!grp.isMassMatrix()) {
    status = grp.computeMassMatrix();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }
  else
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(NOX::Abstract::Group::Ok, 
						   finalStatus,
						   callingFunction);

  status = grp.applyComplex(yVector, zVector, w, *baseCeRealVectorPtr, 
			    *baseCeImagVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Now that Ce is known, call other routine
  status = 
    computeDCeDxa(grp, yVector, zVector, w, aVector, *baseCeRealVectorPtr, 
		  *baseCeImagVectorPtr, result_real, result_imag);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  delete baseCeRealVectorPtr; 
  delete baseCeImagVectorPtr; 

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDCeDxa(
			    LOCA::Bifurcation::HopfBord::AbstractGroup& grp,
			    const NOX::Abstract::Vector& yVector,
			    const NOX::Abstract::Vector& zVector,
			    double w,
			    const NOX::Abstract::MultiVector& aVector,
			    const NOX::Abstract::Vector& Ce_real,
			    const NOX::Abstract::Vector& Ce_imag,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDCeDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Copy original solution vector
  NOX::Abstract::Vector *Xvec = grp.getX().clone(NOX::DeepCopy);

  // Loop over each column of multivector
  for (int i=0; i<aVector.numVectors(); i++) {

    // Perturb solution vector in direction of aVector, return perturbation
    double eps = perturbXVec(grp, *Xvec, aVector[i]);

    // Compute perturbed Ce vectors
    finalStatus = grp.computeJacobian();
    LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

    status = grp.computeMassMatrix();
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

    status = 
      grp.applyComplex(yVector, zVector, w, result_real[i],
		       result_imag[i]);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

    // Difference perturbed and base vector and return approximate derivative
    result_real[i].update(-1.0, Ce_real, 1.0); result_real[i].scale(1.0/eps);
    result_imag[i].update(-1.0, Ce_imag, 1.0); result_imag[i].scale(1.0/eps);

  }

  // Restore original solution vector
  grp.setX(*Xvec);

  delete Xvec;

  return finalStatus;
}

//
// Protected methods start here.
//

double 
LOCA::DerivUtils::perturbParam(LOCA::Continuation::AbstractGroup& grp,
			       double& paramOrig, 
			       int param_id) const
{
  paramOrig = grp.getParam(param_id);

  // Find perturbation size and perturb parameter
  double eps = epsScalar(paramOrig);
  double param = paramOrig + eps;

  // Copy this perturbed parameter vector into group
  grp.setParam(param_id, param);

  // Return perturbation size
  return eps;
}

double 
LOCA::DerivUtils::perturbXVec(LOCA::Continuation::AbstractGroup& grp,
			      const NOX::Abstract::Vector& xVector,
			      const NOX::Abstract::Vector& aVector) const
{
  // Allocate tempertory xVector
  NOX::Abstract::Vector *tmpXVecPtr = xVector.clone(NOX::DeepCopy);

  // Get perturbation size for directional derivative
  double eps = epsVector(*tmpXVecPtr, aVector);

  // Perturb temp vector and copy into group's x vector
  grp.setX(tmpXVecPtr->update(eps, aVector, 1.0));
  delete tmpXVecPtr;

  // Return perturbation size
  return eps;
}

double 
LOCA::DerivUtils::epsScalar(double p) const
{
   return perturb * (perturb + fabs(p));
}

double 
LOCA::DerivUtils::epsVector(const NOX::Abstract::Vector& xVector,
			    const NOX::Abstract::Vector& aVector) const
{
   return perturb * (perturb + xVector.norm(NOX::Abstract::Vector::TwoNorm)
                  / (aVector.norm(NOX::Abstract::Vector::TwoNorm) + perturb));
}
