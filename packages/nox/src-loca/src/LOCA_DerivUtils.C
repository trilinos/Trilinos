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
			      const int param_id, 
			      NOX::Abstract::Vector& result) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDfDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Compute base residual F
  finalStatus = grp.computeF();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  // Allocate new group that we can perturb
  LOCA::Continuation::AbstractGroup* grpPerturbedPtr = 
    dynamic_cast<LOCA::Continuation::AbstractGroup*>(grp.clone());
  
  // Perturb single parameter in this group, and return perturbation, eps
  double eps = perturbParam(*grpPerturbedPtr, param_id);

  // Compute perturbed residual
  status = grpPerturbedPtr->computeF(); 
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Difference perturbed and base vector and return approximate derivative
  NOX::Abstract::Vector *perturbedFPtr = 
    grpPerturbedPtr->getF().clone(NOX::DeepCopy);
  result = doDifference(*perturbedFPtr, grp.getF(), eps);

  delete perturbedFPtr;
  delete grpPerturbedPtr;

  // Compute base residual F
  status = grp.computeF();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  status = grp.computeJacobian();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

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

  finalStatus = grp.computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  status = grp.applyJacobian(nullVector, *baseJnVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Now that Jn is known, call other routine
  status = 
    computeDJnDp(grp, nullVector, param_id, *baseJnVectorPtr, result);
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

  // Form new group that we can perturb
  LOCA::Continuation::AbstractGroup* grpPerturbedPtr = 
    dynamic_cast<LOCA::Continuation::AbstractGroup*>(grp.clone());

  // Perturb single parameter in this group, and return perturbation
  double eps = perturbParam(*grpPerturbedPtr, param_id);

  // Allocate perturbed Jn vector and fill
  NOX::Abstract::Vector *perturbedJnVectorPtr = 
    nullVector.clone(NOX::ShapeCopy);

  finalStatus = grpPerturbedPtr->computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  status = 
    grpPerturbedPtr->applyJacobian(nullVector, *perturbedJnVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Difference perturbed and base vector and return approximate derivative
  result = doDifference(*perturbedJnVectorPtr, JnVector, eps);

  delete perturbedJnVectorPtr; 
  delete grpPerturbedPtr;

  status = grp.computeJacobian();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

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

  // Allocate base Jn vector and fill with J time n
  NOX::Abstract::Vector *baseJnVectorPtr = nullVector.clone(NOX::ShapeCopy);
  
  finalStatus = grp.computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  status = grp.applyJacobian(nullVector, *baseJnVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Now that Jn is known, call other routine
  status = 
    computeDJnDxa(grp, nullVector, aVector, *baseJnVectorPtr, result);
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

  // Form new group that we can perturb
  LOCA::Continuation::AbstractGroup* grpPerturbedPtr = 
    dynamic_cast<LOCA::Continuation::AbstractGroup*>(grp.clone());

  // Perturb solution vector in direction of aVector, return perturbation
  double eps = perturbXVec(*grpPerturbedPtr, aVector);

  // Allocate perturbed Jn vector and fill
  NOX::Abstract::Vector *perturbedJnVectorPtr = 
    nullVector.clone(NOX::ShapeCopy);

  finalStatus = grpPerturbedPtr->computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);
  
  status = 
    grpPerturbedPtr->applyJacobian(nullVector, *perturbedJnVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Difference perturbed and base vector and return approximate derivative
  result = doDifference(*perturbedJnVectorPtr, JnVector, eps);

  delete perturbedJnVectorPtr; 
  delete grpPerturbedPtr;

  status = grp.computeJacobian();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::DerivUtils::computeDJnDxaDp(LOCA::Continuation::AbstractGroup& grp,
				  const NOX::Abstract::Vector& nullVector,
				  const NOX::Abstract::Vector& aVector,
				  const int param_id,
				  NOX::Abstract::Vector& result) const
{
  string callingFunction = 
    "LOCA::DerivUtils::computeDJnDxaDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Allocate base Jn vector and fill with J time n
  NOX::Abstract::Vector *baseJnVectorPtr = nullVector.clone(NOX::ShapeCopy);
  
  finalStatus = grp.computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  status = grp.applyJacobian(nullVector, *baseJnVectorPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Now that Jn is known, call other routine
  status = 
    computeDJnDxa(grp, nullVector, aVector, *baseJnVectorPtr, result);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  NOX::Abstract::Vector *resultVecPtr = nullVector.clone(NOX::ShapeCopy);
  status = computeDJnDp(grp, nullVector, param_id, *baseJnVectorPtr, 
			*resultVecPtr);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  result.update(1.0, *resultVecPtr, 1.0);

  delete baseJnVectorPtr; 
  delete resultVecPtr; 

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

  finalStatus = grp.computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  status = grp.computeMassMatrix();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
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

  // Form new group that we can perturb
  LOCA::Bifurcation::HopfBord::AbstractGroup* grpPerturbedPtr = 
    dynamic_cast<LOCA::Bifurcation::HopfBord::AbstractGroup*>(grp.clone());

  // Perturb single parameter in this group, and return perturbation
  double eps = perturbParam(*grpPerturbedPtr, param_id);

  // Compute perturbed Ce vectors
  finalStatus = grpPerturbedPtr->computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  status = grpPerturbedPtr->computeMassMatrix();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  status = 
    grpPerturbedPtr->applyComplex(yVector, zVector, w, result_real,
				  result_imag);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Difference perturbed and base vector and return approximate derivative
  doDifference2(result_real, Ce_real, eps);
  doDifference2(result_imag, Ce_imag, eps);

  delete grpPerturbedPtr;

  grp.computeJacobian();
  grp.computeMassMatrix();

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

  finalStatus = grp.computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  status = grp.computeMassMatrix();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
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

  // Form new group that we can perturb
  LOCA::Bifurcation::HopfBord::AbstractGroup* grpPerturbedPtr = 
    dynamic_cast<LOCA::Bifurcation::HopfBord::AbstractGroup*>(grp.clone());

  // Perturb solution vector in direction of aVector, return perturbation
  double eps = perturbXVec(*grpPerturbedPtr, aVector);

  // Compute perturbed Ce vectors
  finalStatus = grpPerturbedPtr->computeJacobian();
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  status = grpPerturbedPtr->computeMassMatrix();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  status = 
    grpPerturbedPtr->applyComplex(yVector, zVector, w, result_real,
				  result_imag);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Difference perturbed and base vector and return approximate derivative
  doDifference2(result_real, Ce_real, eps);
  doDifference2(result_imag, Ce_imag, eps);

  delete grpPerturbedPtr;

  status = grp.computeJacobian();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  status = grp.computeMassMatrix();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  return finalStatus;
}

//
// Protected methods start here.
//

NOX::Abstract::Vector& 
LOCA::DerivUtils::doDifference(NOX::Abstract::Vector& perturbedVector,
			       const NOX::Abstract::Vector& baseVector,
			       const double eps) const
{
  perturbedVector.update(-1.0, baseVector, 1.0);

  return perturbedVector.scale(1.0/eps);
}

void
LOCA::DerivUtils::doDifference2(NOX::Abstract::Vector& perturbedVector,
			       const NOX::Abstract::Vector& baseVector,
			       const double eps) const
{
  perturbedVector.update(-1.0, baseVector, 1.0);
  perturbedVector.scale(1.0/eps);
}

double 
LOCA::DerivUtils::perturbParam(
			LOCA::Continuation::AbstractGroup& grpPerturbed,
			const int param_id) const
{
  // Allocated new parameter vector
  LOCA::ParameterVector paramVec(grpPerturbed.getParams());

  // Find perturbation size and perturb parameter
  double eps = epsScalar(paramVec[param_id]);
  paramVec[param_id] += eps;

  // Copy this perturbed parameter vector into group
  grpPerturbed.setParams(paramVec);

  // Return perturbation size
  return eps;
}

double 
LOCA::DerivUtils::perturbXVec(LOCA::Continuation::AbstractGroup& grpPerturbed,
			      const NOX::Abstract::Vector& aVector) const
{
  // Allocate tempertory xVector
  NOX::Abstract::Vector *tmpXVecPtr = grpPerturbed.getX().clone(NOX::DeepCopy);

  // Get perturbation size for directional derivative
  double eps = epsVector(*tmpXVecPtr, aVector);

  // Perturb temp vector and copy into group's x vector
  grpPerturbed.setX(tmpXVecPtr->update(eps, aVector, 1.0));
  delete tmpXVecPtr;

  // Return perturbation size
  return eps;
}

double 
LOCA::DerivUtils::epsScalar(const double p) const
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
