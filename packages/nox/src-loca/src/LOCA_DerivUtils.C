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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER


// Class Header 
#include "LOCA_DerivUtils.H"

#include "LOCA_Abstract_Group.H" 
#include "NOX_Abstract_Vector.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Common.H"  // For fabs function

using namespace LOCA;

DerivUtils::DerivUtils(double perturb) :
	perturb(perturb)
{
  // assert (perturb>0.0);
}

DerivUtils::DerivUtils(const DerivUtils& source) :
	perturb(source.perturb)
{

}

DerivUtils::~DerivUtils()
{
 
}

DerivUtils* DerivUtils::clone(NOX::CopyType type)
{
  return new DerivUtils(*this);  //Call Copy Constructor
}

bool DerivUtils::computeDfDp(LOCA::Abstract::Group& grp,
       const int param_id, NOX::Abstract::Vector& result) const
{
  // Compute base residual F
  grp.computeF();

  // Allocate new group that we can perturb
  LOCA::Abstract::Group* grpPerturbedPtr = grp.clone();
  
  // Perturb single parameter in this group, and return perturbation, eps
  double eps = perturbParam(*grpPerturbedPtr, param_id);

  // Compute perturbed residual
  bool ok = grpPerturbedPtr->computeF(); 

  // Difference perturbed and base vector and return approximate derivative
  NOX::Abstract::Vector *perturbedFPtr = grpPerturbedPtr->getF().clone(NOX::DeepCopy);
  result = doDifference(*perturbedFPtr, grp.getF(), eps);

  delete perturbedFPtr;
  delete grpPerturbedPtr;

  return ok;
}

bool DerivUtils::computeDJnDp(LOCA::Abstract::Group& grp,
       const NOX::Abstract::Vector& nullVector, const int param_id,
                               NOX::Abstract::Vector& result) const
{
  // Allocate base Jn vector and fill with J time n
  NOX::Abstract::Vector *baseJnVectorPtr = nullVector.clone(NOX::ShapeCopy);
  grp.computeJacobian();
  grp.applyJacobian(nullVector, *baseJnVectorPtr);

  // Now that Jn is known, call other routine
  bool ok = computeDJnDp(grp, nullVector, param_id, *baseJnVectorPtr, result);

  delete baseJnVectorPtr; 

  return ok;
}

bool DerivUtils::computeDJnDp(LOCA::Abstract::Group& grp,
       const NOX::Abstract::Vector& nullVector, const int param_id,
       const NOX::Abstract::Vector& JnVector,
       NOX::Abstract::Vector& result) const
{
  // Form new group that we can perturb
  LOCA::Abstract::Group* grpPerturbedPtr = grp.clone();

  // Perturb single parameter in this group, and return perturbation
  double eps = perturbParam(*grpPerturbedPtr, param_id);

  // Allocate perturbed Jn vector and fill
  NOX::Abstract::Vector *perturbedJnVectorPtr = nullVector.clone(NOX::ShapeCopy);
  grpPerturbedPtr->computeJacobian();
  bool ok = grpPerturbedPtr->applyJacobian(nullVector, *perturbedJnVectorPtr);

  // Difference perturbed and base vector and return approximate derivative
  result = doDifference(*perturbedJnVectorPtr, JnVector, eps);

  delete perturbedJnVectorPtr; 
  delete grpPerturbedPtr;

  return ok;
}

bool DerivUtils::computeDJnDxa(LOCA::Abstract::Group& grp,
                const NOX::Abstract::Vector& nullVector,
                const NOX::Abstract::Vector& aVector,
                NOX::Abstract::Vector& result) const
{
  // Allocate base Jn vector and fill with J time n
  NOX::Abstract::Vector *baseJnVectorPtr = nullVector.clone(NOX::ShapeCopy);
  grp.computeJacobian();
  grp.applyJacobian(nullVector, *baseJnVectorPtr);

  // Now that Jn is known, call other routine
  bool ok = computeDJnDxa(grp, nullVector, aVector, *baseJnVectorPtr, result);

  delete baseJnVectorPtr; 

  return ok;
}

bool DerivUtils::computeDJnDxa(LOCA::Abstract::Group& grp,
                 const NOX::Abstract::Vector& nullVector,
                 const NOX::Abstract::Vector& aVector,
                 const NOX::Abstract::Vector& JnVector,
                 NOX::Abstract::Vector& result) const
{
  // Form new group that we can perturb
  LOCA::Abstract::Group* grpPerturbedPtr = grp.clone();

  // Perturb solution vector in direction of aVector, return perturbation
  double eps = perturbXVec(*grpPerturbedPtr, aVector);

  // Allocate perturbed Jn vector and fill
  NOX::Abstract::Vector *perturbedJnVectorPtr = nullVector.clone(NOX::ShapeCopy);
  grpPerturbedPtr->computeJacobian();
  bool ok = grpPerturbedPtr->applyJacobian(nullVector, *perturbedJnVectorPtr);

  // Difference perturbed and base vector and return approximate derivative
  result = doDifference(*perturbedJnVectorPtr, JnVector, eps);

  delete perturbedJnVectorPtr; 
  delete grpPerturbedPtr;

  return ok;
}

bool DerivUtils::computeDJnDxaDp(LOCA::Abstract::Group& grp,
                             const NOX::Abstract::Vector& nullVector,
                             const NOX::Abstract::Vector& aVector,
                             const int param_id,
                             NOX::Abstract::Vector& result) const
{
  // Allocate base Jn vector and fill with J time n
  NOX::Abstract::Vector *baseJnVectorPtr = nullVector.clone(NOX::ShapeCopy);
  grp.computeJacobian();
  grp.applyJacobian(nullVector, *baseJnVectorPtr);

  // Now that Jn is known, call other routine
  bool ok = computeDJnDxa(grp, nullVector, aVector, *baseJnVectorPtr, result);

  NOX::Abstract::Vector *resultVecPtr = nullVector.clone(NOX::ShapeCopy);
  computeDJnDp(grp, nullVector, param_id, *baseJnVectorPtr, *resultVecPtr);

  result.update(1.0, *resultVecPtr, 1.0);

  delete baseJnVectorPtr; 
  delete resultVecPtr; 

  return ok;
}

//
// Protected methods start here.
//

NOX::Abstract::Vector& DerivUtils::doDifference(
                     NOX::Abstract::Vector& perturbedVector,
                     const NOX::Abstract::Vector& baseVector,
                     const double eps) const
{
  perturbedVector.update(-1.0, baseVector, 1.0);

  return perturbedVector.scale(1.0/eps);
}

double DerivUtils::perturbParam(LOCA::Abstract::Group& grpPerturbed,
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

double DerivUtils::perturbXVec(LOCA::Abstract::Group& grpPerturbed,
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

double DerivUtils::epsScalar(const double p) const
{
   return perturb * (perturb + fabs(p));
}

double DerivUtils::epsVector(const NOX::Abstract::Vector& xVector,
                                const NOX::Abstract::Vector& aVector) const
{
   return perturb * (perturb + xVector.norm(NOX::Abstract::Vector::TwoNorm)
                               / aVector.norm(NOX::Abstract::Vector::TwoNorm));
}
