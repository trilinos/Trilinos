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

#include "LOCA_Bifurcation_ArcLengthGroup.H"
#include "LOCA_Parameter_Vector.H"

using namespace LOCA;
using namespace Bifurcation;
using NOX::Abstract::Vector;

LOCA::Bifurcation::ArcLengthGroup::ArcLengthGroup(const Abstract::Group& g,
					    int paramId,
					    double param,
                                            const DerivUtils& d)
  : grpPtr(dynamic_cast<LOCA::Abstract::Group *>(g.clone(NOX::DeepCopy))), 
    alXVec(g.getX(), param),
    alFVec(g.getX(), 0.0),
    alNewtonVec(g.getX(), 0.0),
    alPrevXVec(g.getX(), param),
    alTangentVec(g.getX(), 0.0),
    derivResidualParamPtr(g.getX().clone(NOX::ShapeCopy)),
    derivPtr(d.clone()),
    arcParamId(paramId), 
    arclengthStep(0.0),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false) {}

LOCA::Bifurcation::ArcLengthGroup::ArcLengthGroup(const ArcLengthGroup& source)
  : grpPtr(dynamic_cast<LOCA::Abstract::Group *>
          (source.grpPtr->clone(NOX::DeepCopy))), 
    alXVec(source.alXVec),
    alFVec(source.alFVec),
    alNewtonVec(source.alNewtonVec),
    alPrevXVec(source.alPrevXVec),
    alTangentVec(source.alTangentVec),
    derivResidualParamPtr(source.grpPtr->getX().clone(NOX::ShapeCopy)),
    derivPtr(source.derivPtr->clone()),
    arcParamId(source.arcParamId),
    arclengthStep(source.arclengthStep),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton) {}


LOCA::Bifurcation::ArcLengthGroup::~ArcLengthGroup() 
{
  delete derivResidualParamPtr;
  delete derivPtr;
  delete grpPtr;
}

LOCA::Abstract::Group&
LOCA::Bifurcation::ArcLengthGroup::operator=(const LOCA::Abstract::Group& source)
{
  return *this = dynamic_cast<const ArcLengthGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Bifurcation::ArcLengthGroup::operator=(const NOX::Abstract::Group& source)
{
  return *this = dynamic_cast<const ArcLengthGroup&>(source);
}

ArcLengthGroup&
LOCA::Bifurcation::ArcLengthGroup::operator=(const ArcLengthGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    NOX::CopyType type = NOX::DeepCopy;

    // Delete old values
    delete grpPtr;

    // Copy values
    grpPtr = dynamic_cast<LOCA::Abstract::Group *>(source.grpPtr->clone(type));
    alXVec = source.alXVec;
    alFVec = source.alFVec;
    alNewtonVec = source.alNewtonVec;
    alPrevXVec = source.alPrevXVec;
    alTangentVec = source.alTangentVec;
    arcParamId = source.arcParamId;
    arclengthStep = source.arclengthStep;
    derivResidualParamPtr= source.derivResidualParamPtr;
    derivPtr = source.derivPtr;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
  }

  return *this;
}

NOX::Abstract::Group*
LOCA::Bifurcation::ArcLengthGroup::clone(NOX::CopyType type) const 
{
  return new ArcLengthGroup(*this);
}

void
LOCA::Bifurcation::ArcLengthGroup::setParams(const ParameterVector& p) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(p);
}

void
LOCA::Bifurcation::ArcLengthGroup::computeParams(const ParameterVector& oldParams,
					      const ParameterVector& direction,
					      double step) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->computeParams(oldParams, direction, step);
}

const ParameterVector&
LOCA::Bifurcation::ArcLengthGroup::getParams() const 
{
  return grpPtr->getParams();
}

void
LOCA::Bifurcation::ArcLengthGroup::setArclengthStep(double deltaS) 
{
  arclengthStep = deltaS;
  
  // Arclength Step appears on RHS but not in Jacobian
  isValidF = false;
  isValidNewton = false;
}

double
LOCA::Bifurcation::ArcLengthGroup::getArcParam() const 
{
  ParameterVector params(grpPtr->getParams());
  return params[arcParamId];
}

void
LOCA::Bifurcation::ArcLengthGroup::setArcParam(double param) 
{
  ParameterVector params(grpPtr->getParams());

  params[arcParamId] = param;
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(params);
}

void 
LOCA::Bifurcation::ArcLengthGroup::setX(const Vector& y) 
{
  setX( dynamic_cast<const ArcLengthVector&>(y) );
}

void
LOCA::Bifurcation::ArcLengthGroup::setX(const ArcLengthVector& y) 
{
  grpPtr->setX( y.getXVec() );
  alXVec = y;

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Bifurcation::ArcLengthGroup::setPrevX(const Vector& y) 
{
  setPrevX( dynamic_cast<const ArcLengthVector&>(y) );
}

void
LOCA::Bifurcation::ArcLengthGroup::setPrevX(const ArcLengthVector& y) 
{
  alPrevXVec = y;

  isValidF = false; // Previous vector part of residual equation
  isValidNewton = false;
}

void
LOCA::Bifurcation::ArcLengthGroup::computeX(const NOX::Abstract::Group& g, 
					 const Vector& d,
					 double step) 
{
  computeX( dynamic_cast<const ArcLengthGroup&>(g),
	    dynamic_cast<const ArcLengthVector&>(d),
	    step);
}

void
LOCA::Bifurcation::ArcLengthGroup::computeX(const ArcLengthGroup& g, 
					    const ArcLengthVector& d,
					    double step) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->computeX(*(g.grpPtr), d.getXVec(), step);
  alXVec.update(1.0, g.getX(), step, d, 0.0);
  setArcParam(alXVec.getArcParam()); 
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::ArcLengthGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res;

  res = grpPtr->computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  alFVec.getXVec() = grpPtr->getF();
  
  // Construct residual of arclength equation
  ArcLengthVector *tmpVec =
    dynamic_cast<ArcLengthVector *>(alPrevXVec.clone(NOX::DeepCopy));
  tmpVec->update(1.0, alXVec, -1.0);
  
  alFVec.getArcParam() = alTangentVec.dot(*tmpVec) - arclengthStep;

  delete tmpVec;
  
  isValidF = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::ArcLengthGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = grpPtr->computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  res = derivPtr->computeDfDp(*grpPtr, arcParamId, 
			      *derivResidualParamPtr);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  isValidJacobian = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::ArcLengthGroup::computeGradient() 
{
  return NOX::Abstract::Group::NotDefined;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::ArcLengthGroup::computeNewton(NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  res = computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = applyJacobianInverse(params, alFVec, alNewtonVec);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  alNewtonVec.scale(-1.0);

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::ArcLengthGroup::computeTangent(NOX::Parameter::List& params,
                                                  int paramId) 
{
  Vector& tangent_x = alTangentVec.getXVec();
  double& tangent_param =  alTangentVec.getArcParam();

  Vector *dfdpVec = tangent_x.clone(NOX::ShapeCopy);

  NOX::Abstract::Group::ReturnType res = grpPtr->computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = grpPtr->computeDfDp(arcParamId, *dfdpVec);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = grpPtr->applyJacobianInverse(params, *dfdpVec, tangent_x);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  tangent_param = sqrt(1.0 + tangent_x.dot(tangent_x));
  
  // Construct residual of arclength equation
  ArcLengthVector tmpVec(alPrevXVec);
  tmpVec.update(1.0, alXVec, -1.0);
   
  ArcLengthVector tmpVec2(*dfdpVec, 1.0);

  // CHECK FOR SIGN CHANGE!
  if (tmpVec.dot(tmpVec2) < 0.0) tangent_param *= -1.0;

  tangent_x.scale(tangent_param);

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  alTangentVec.setVec(tangent_x, tangent_param);

  delete dfdpVec;
  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::ArcLengthGroup::applyJacobian(const Vector& input,
					      Vector& result) const 
{
  const ArcLengthVector& al_input = dynamic_cast<const ArcLengthVector&>(input);
  ArcLengthVector& al_result = dynamic_cast<ArcLengthVector&>(result);

  const Vector& input_x = al_input.getXVec();
  double input_param = al_input.getArcParam();

  Vector *result_x = input_x.clone(NOX::ShapeCopy);
  double result_param;

  double arcParam = al_input.getArcParam();  //ADDED al_input. 

  NOX::Abstract::Group::ReturnType res;

  // compute J*x
  res = grpPtr->applyJacobian(input_x, *result_x);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // compute J*x + p*dR/dp
  result_x->update(arcParam, *derivResidualParamPtr, 1.0);

  // compute dx/ds x + dp/ds p
  result_param = alTangentVec.dot(al_input);

  al_result.setVec(*result_x, result_param);

  delete result_x;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::ArcLengthGroup::applyJacobianTranspose(const Vector& input,
						       Vector& result) const 
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::ArcLengthGroup::applyJacobianInverse(NOX::Parameter::List& params,
						     const Vector& input,
						     Vector& result) const 
{
  const ArcLengthVector& al_input = dynamic_cast<const ArcLengthVector&>(input);
  ArcLengthVector& al_result = dynamic_cast<ArcLengthVector&>(result);

  const Vector& input_x = al_input.getXVec();
  double input_param = al_input.getArcParam();

  Vector *result_x = input_x.clone(NOX::ShapeCopy);
  double result_param;

  Vector *a = input_x.clone(NOX::ShapeCopy);
  Vector *b = input_x.clone(NOX::ShapeCopy);
  
  NOX::Abstract::Group::ReturnType res;

  // Solve J*a = input_x
  if (input_x.norm() > 0.0) {
    res = grpPtr->applyJacobianInverse(params, input_x, *a);
    if (res != NOX::Abstract::Group::Ok)
      return res;
  }

  // Solve J*b = dR/dp
  res = grpPtr->applyJacobianInverse(params, *derivResidualParamPtr, *b);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute result_param = 
  result_param =  -(input_param - alTangentVec.getXVec().dot(*a))
                / (alTangentVec.getXVec().dot(*b) - alTangentVec.getArcParam());

  // Compute result_x = a + result_param*b 
  *result_x = result_x->update(1.0, *a, -result_param, *b, 0.0);

  al_result.setVec(*result_x, result_param);

  // Clean up memory
  delete result_x;
  delete a;
  delete b;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::ArcLengthGroup::applyJacobianDiagonalInverse(const Vector& input, Vector& result) const 
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::ArcLengthGroup::applyRightPreconditioning(NOX::Parameter::List& params, const Vector& input, Vector& result) const 
{
  const ArcLengthVector& al_input = dynamic_cast<const ArcLengthVector&>(input);
  ArcLengthVector& al_result = dynamic_cast<ArcLengthVector&>(result);

  const Vector& input_x = al_input.getXVec();
  double input_param = al_input.getArcParam();

  Vector *result_x = input_x.clone(NOX::ShapeCopy);
  double result_param;

  Vector *a = input_x.clone(NOX::ShapeCopy);
  Vector *b = input_x.clone(NOX::ShapeCopy);

  NOX::Abstract::Group::ReturnType res;

  // Solve J*a = input_x
  res = grpPtr->applyRightPreconditioning(params, input_x, *a);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Solve J*b = dR/dp
  res = grpPtr->applyRightPreconditioning(params, *derivResidualParamPtr, *b);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute result_param = 
  result_param =  (input_param - alTangentVec.getXVec().dot(*a))
                / (alTangentVec.getXVec().dot(*b) - alTangentVec.getArcParam());

  // Compute result_x = a + result_param*b 
  *result_x = result_x->update(1.0, *a, result_param, *b, 0.0);

  al_result.setVec(*result_x, result_param);

  // Clean up memory
  delete result_x;
  delete a;
  delete b;

  return res;
}

bool
LOCA::Bifurcation::ArcLengthGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::Bifurcation::ArcLengthGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Bifurcation::ArcLengthGroup::isGradient() const 
{
  return false;
}

bool
LOCA::Bifurcation::ArcLengthGroup::isNewton() const 
{
  return isValidNewton;
}
  
const Vector&
LOCA::Bifurcation::ArcLengthGroup::getX() const 
{
  return alXVec;
}

const Vector&
LOCA::Bifurcation::ArcLengthGroup::getF() const 
{
  return alFVec;
}

double
LOCA::Bifurcation::ArcLengthGroup::getNormF() const 
{
  return alFVec.norm();
}

const Vector&
LOCA::Bifurcation::ArcLengthGroup::getGradient() const 
{
  return getNewton();
}

const Vector&
LOCA::Bifurcation::ArcLengthGroup::getNewton() const 
{
  return alNewtonVec;
}

const LOCA::Bifurcation::ArcLengthVector&
LOCA::Bifurcation::ArcLengthGroup::getTangent() const 
{
  return alTangentVec;
}

double
LOCA::Bifurcation::ArcLengthGroup::getNormNewtonSolveResidual() const 
{
  ArcLengthVector residual = alFVec;
  
  NOX::Abstract::Group::ReturnType res = applyJacobian(alNewtonVec, residual);
  if (res != NOX::Abstract::Group::Ok) {
    cout << "ERROR: applyJacobian() in getNormNewtonSolveResidual "
	 << " returned not ok" << endl;
    throw "LOCA Error";
    return 0.0;
  }

  residual = residual.update(1.0, alFVec, 1.0);
  return residual.norm();
}

void
LOCA::Bifurcation::ArcLengthGroup::print() const
{
  cout << "Beginning ArcLengthGroup.print:" << endl;
  cout << "Underlying Group = " << endl;
  grpPtr->print();

  cout << endl;

  cout << "alXVec = " << endl;
  alXVec.print();

  cout << endl;

  if (isValidF) {
    cout << "alFVec = " << endl;
    alFVec.print();
  }
  else
    cout << "alFVec not computed" << endl;

  cout << endl;
}

