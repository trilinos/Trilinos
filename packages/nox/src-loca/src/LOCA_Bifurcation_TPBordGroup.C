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

#include "LOCA_Bifurcation_TPBordGroup.H"
#include "LOCA_Parameter_Vector.H"

using namespace LOCA;
using namespace Bifurcation;
using NOX::Abstract::Vector;

LOCA::Bifurcation::TPBordGroup::TPBordGroup(const Abstract::Group& g,
					    const Vector& lenVec,
					    int paramId,
					    const DerivUtils& d)
  : grpPtr(dynamic_cast<Abstract::Group*>(g.clone(NOX::DeepCopy))), 
    tpXVec(g.getX(), lenVec, 0.0),
    tpFVec(lenVec, lenVec, 0.0),
    tpNewtonVec(lenVec, lenVec, 0.0),
    tpTangentVec(lenVec, lenVec, 0.0),
    lengthVecPtr(lenVec.clone(NOX::ShapeCopy)), 
    bifParamId(paramId), 
    derivResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    derivNullResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    derivPtr(d.clone(NOX::DeepCopy)), 
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false) 
{
  tpXVec.getBifParam() = getBifParam();
}

LOCA::Bifurcation::TPBordGroup::TPBordGroup(const TPBordGroup& source)
  : grpPtr(dynamic_cast<Abstract::Group*>(source.grpPtr->clone(NOX::DeepCopy))), 
    tpXVec(source.tpXVec),
    tpFVec(source.tpFVec),
    tpNewtonVec(source.tpNewtonVec),
    tpTangentVec(source.tpTangentVec),
    lengthVecPtr(source.lengthVecPtr->clone(NOX::DeepCopy)), 
    bifParamId(source.bifParamId),
    derivResidualParamPtr(source.derivResidualParamPtr->clone(NOX::DeepCopy)),
    derivNullResidualParamPtr(source.derivNullResidualParamPtr->clone(NOX::DeepCopy)),
    derivPtr(source.derivPtr->clone(NOX::DeepCopy)),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton) {}


LOCA::Bifurcation::TPBordGroup::~TPBordGroup() 
{
  delete grpPtr;
  delete lengthVecPtr;
  delete derivResidualParamPtr;
  delete derivNullResidualParamPtr;
  delete derivPtr;
  
}

LOCA::Abstract::Group&
LOCA::Bifurcation::TPBordGroup::operator=(const LOCA::Abstract::Group& source)
{
  return *this = dynamic_cast<const TPBordGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Bifurcation::TPBordGroup::operator=(const NOX::Abstract::Group& source)
{
  return *this = dynamic_cast<const TPBordGroup&>(source);
}

TPBordGroup&
LOCA::Bifurcation::TPBordGroup::operator=(const TPBordGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    NOX::CopyType type = NOX::DeepCopy;

    // Delete old values
    delete grpPtr;
    delete lengthVecPtr;
    delete derivResidualParamPtr;
    delete derivNullResidualParamPtr;
    delete derivPtr;

    // Copy values
    grpPtr = dynamic_cast<Abstract::Group*>(source.grpPtr->clone(type));
    tpXVec = source.tpXVec;
    tpFVec = source.tpFVec;
    tpNewtonVec = source.tpNewtonVec;
    tpTangentVec = source.tpTangentVec;
    lengthVecPtr = source.lengthVecPtr->clone(type);
    derivResidualParamPtr = source.derivResidualParamPtr->clone(type);
    derivNullResidualParamPtr = source.derivNullResidualParamPtr->clone(type);
    derivPtr = source.derivPtr->clone(type);
    bifParamId = source.bifParamId;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
  }

  return *this;
}

NOX::Abstract::Group*
LOCA::Bifurcation::TPBordGroup::clone(NOX::CopyType type) const 
{
  return new TPBordGroup(*this);
}

bool
LOCA::Bifurcation::TPBordGroup::setParams(const ParameterVector& p) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  return grpPtr->setParams(p);
}

bool
LOCA::Bifurcation::TPBordGroup::computeParams(const ParameterVector& oldParams,
					      const ParameterVector& direction,
					      double step) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  return grpPtr->computeParams(oldParams, direction, step);
}

const ParameterVector&
LOCA::Bifurcation::TPBordGroup::getParams() const 
{
  return grpPtr->getParams();
}

double
LOCA::Bifurcation::TPBordGroup::getBifParam() const 
{
  ParameterVector params(grpPtr->getParams());
  return params[bifParamId];
}

bool
LOCA::Bifurcation::TPBordGroup::setBifParam(double param) 
{
  ParameterVector params(grpPtr->getParams());

  params[bifParamId] = param;
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  return grpPtr->setParams(params);
}

bool
LOCA::Bifurcation::TPBordGroup::setX(const Vector& y) 
{
  return setX( dynamic_cast<const TPBordVector&>(y) );
}

bool
LOCA::Bifurcation::TPBordGroup::setX(const TPBordVector& y) 
{
  bool res;

  res = grpPtr->setX( y.getXVec() );
  tpXVec = y;

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  return res;
}

bool
LOCA::Bifurcation::TPBordGroup::computeX(const NOX::Abstract::Group& g, 
					 const Vector& d,
					 double step) 
{
  return computeX( dynamic_cast<const TPBordGroup&>(g),
		   dynamic_cast<const TPBordVector&>(d),
		   step);
}

bool
LOCA::Bifurcation::TPBordGroup::computeX(const TPBordGroup& g, 
					 const TPBordVector& d,
					 double step) 
{
  bool res;

  res = grpPtr->computeX(*(g.grpPtr), d.getXVec(), step);
  tpXVec.update(1.0, g.getX(), step, d, 0.0);
  res = res & setBifParam(tpXVec.getBifParam());

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  return res;
}

bool
LOCA::Bifurcation::TPBordGroup::computeF() 
{
  bool res = true;

  if (isValidF == false) {
    res = res & grpPtr->computeF();
    tpFVec.getXVec() = grpPtr->getF();

    res = res & grpPtr->computeJacobian();
    res = res & grpPtr->applyJacobian(tpXVec.getNullVec(), 
				      tpFVec.getNullVec());

    tpFVec.getBifParam() = tpXVec.getNullVec().dot(*lengthVecPtr) - 1.0;
    
    isValidF = true;
    
  }

  return res;
}

bool
LOCA::Bifurcation::TPBordGroup::computeJacobian() 
{
  bool res = computeF();

  if (isValidJacobian == false) {
    res = grpPtr->computeJacobian();

    res = res & derivPtr->computeDfDp(*grpPtr, bifParamId, 
				      *derivResidualParamPtr);
    res = res & derivPtr->computeDJnDp(*grpPtr, tpXVec.getNullVec(), 
				       bifParamId,
				       tpFVec.getNullVec(), 
				       *derivNullResidualParamPtr);
    isValidJacobian = true;
  }

  return res;
}

bool
LOCA::Bifurcation::TPBordGroup::computeGradient() 
{
  return false;
}
   
bool
LOCA::Bifurcation::TPBordGroup::computeNewton(NOX::Parameter::List& params) 
{
  bool res = computeF();
  
  res = res & computeJacobian();

  if (isValidNewton == false) {
    applyJacobianInverse(params, tpFVec, tpNewtonVec);
    tpNewtonVec.scale(-1.0);
    isValidNewton = true;
  }

  return res;
}

bool
LOCA::Bifurcation::TPBordGroup::computeTangent(NOX::Parameter::List& params,
                                      int paramID)
{
  bool res = computeJacobian();

  NOX::Abstract::Vector* dfdpVec = tpTangentVec.clone(NOX::ShapeCopy);

  res = res && computeDfDp(paramID, *dfdpVec);

  res = res && applyJacobianInverse(params, *dfdpVec, tpTangentVec);

  delete dfdpVec;

  return res;
}

bool
LOCA::Bifurcation::TPBordGroup::applyJacobian(const Vector& input,
					      Vector& result) const 
{
  const TPBordVector& tp_input = dynamic_cast<const TPBordVector&>(input);
  TPBordVector& tp_result = dynamic_cast<TPBordVector&>(result);

  const Vector& input_x = tp_input.getXVec();
  const Vector& input_null = tp_input.getNullVec();
  double input_param = tp_input.getBifParam();

  Vector *result_x = input_x.clone(NOX::ShapeCopy);
  Vector *result_null = input_null.clone(NOX::ShapeCopy);
  Vector *tmp = input_null.clone(NOX::ShapeCopy);
  double result_param;

  double bifParam = getBifParam();

  bool res;

  // compute J*x
  res = grpPtr->applyJacobian(input_x, *result_x);

  // compute J*x + p*dR/dp
  result_x->update(bifParam, *derivResidualParamPtr, 1.0);

  // compute J*y
  res = res & grpPtr->applyJacobian(input_null, *result_null);

  // compute J*y + p*dJy/dp
  result_null->update(bifParam, *derivNullResidualParamPtr, 1.0);

  // compute (dJy/dx)*x
  res = res & derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), input_x, 
				      tpFVec.getNullVec(), *tmp);
  
  // compute (dJy/dx)*x + J*y + p*dJy/dp
  result_null->update(1.0, *tmp, 1.0);

  // compute phi^T*y
  result_param = lengthVecPtr->dot(input_null);

  tp_result.setVec(*result_x, *result_null, result_param);

  delete result_x;
  delete result_null;
  delete tmp;

  return res;
}

bool
LOCA::Bifurcation::TPBordGroup::applyJacobianTranspose(const Vector& input,
						       Vector& result) const 
{
  return false;
}

bool
LOCA::Bifurcation::TPBordGroup::applyJacobianInverse(NOX::Parameter::List& params,
						     const Vector& input,
						     Vector& result) const 
{
  const TPBordVector& tp_input = dynamic_cast<const TPBordVector&>(input);
  TPBordVector& tp_result = dynamic_cast<TPBordVector&>(result);

  const Vector& input_x = tp_input.getXVec();
  const Vector& input_null = tp_input.getNullVec();
  double input_param = tp_input.getBifParam();

  Vector *result_x = input_x.clone(NOX::ShapeCopy);
  Vector *result_null = input_null.clone(NOX::ShapeCopy);
  Vector *tmp = input_null.clone(NOX::ShapeCopy);
  double result_param;

  Vector *a = input_x.clone(NOX::ShapeCopy);
  Vector *b = input_x.clone(NOX::ShapeCopy);
  Vector *c = input_x.clone(NOX::ShapeCopy);
  Vector *d = input_x.clone(NOX::ShapeCopy);
  Vector *derivJ = input_x.clone(NOX::ShapeCopy);
  
  bool res;

  // Solve J*a = input_x
  res = res & grpPtr->applyJacobianInverse(params, input_x, *a);

  // Solve J*b = dR/dp
  res = res & grpPtr->applyJacobianInverse(params, *derivResidualParamPtr, *b);

  // Compute (dJy/dx)*a
  res = res & derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), *a, 
				      tpFVec.getNullVec(), *derivJ);

  // Compute input_null - (dJy/dx)*a
  *derivJ = derivJ->update(1.0, input_null, -1.0);

  // Solve J*c = input_null - (dJy/dx)*a
  res = res & grpPtr->applyJacobianInverse(params, *derivJ, *c);

  // Compute (dJy/dx)*b
  res = res & derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), *b, 
				      tpFVec.getNullVec(), *derivJ);

  // Compute (dJy/dx)*b - dJy/dp
  *derivJ = derivJ->update(-1.0, *derivNullResidualParamPtr, 1.0);

  // Solve J*d = (dJy/dx)*b - dJy/dp
  res = res & grpPtr->applyJacobianInverse(params, *derivJ, *d);

  // Compute result_param = -(input_param - phi^T*c)/(phi^T*d)
  result_param = (input_param- lengthVecPtr->dot(*c))/(lengthVecPtr->dot(*d));

  // Compute result_x = a - result_param*b 
  *result_x = result_x->update(1.0, *a, -result_param, *b, 0.0);

  // Compute result_null = c + result_param*d
  *result_null = result_null->update(1.0, *c, result_param, *d, 0.0);

  tp_result.setVec(*result_x, *result_null, result_param);

  // Clean up memory
  delete result_x;
  delete result_null;
  delete tmp;

  delete a;
  delete b;
  delete c;
  delete d;
  delete derivJ;

  return res;

}

bool
LOCA::Bifurcation::TPBordGroup::applyJacobianDiagonalInverse(const Vector& input, Vector& result) const 
{
  return false;
}

bool
LOCA::Bifurcation::TPBordGroup::applyRightPreconditioning(NOX::Parameter::List& params, const Vector& input, Vector& result) const 
{
  const TPBordVector& tp_input = dynamic_cast<const TPBordVector&>(input);
  TPBordVector& tp_result = dynamic_cast<TPBordVector&>(result);

  const Vector& input_x = tp_input.getXVec();
  const Vector& input_null = tp_input.getNullVec();
  double input_param = tp_input.getBifParam();

  Vector *result_x = input_x.clone(NOX::ShapeCopy);
  Vector *result_null = input_null.clone(NOX::ShapeCopy);
  Vector *tmp = input_null.clone(NOX::ShapeCopy);
  double result_param;

  Vector *a = input_x.clone(NOX::ShapeCopy);
  Vector *b = input_x.clone(NOX::ShapeCopy);
  Vector *c = input_x.clone(NOX::ShapeCopy);
  Vector *d = input_x.clone(NOX::ShapeCopy);
  Vector *derivJ = input_x.clone(NOX::ShapeCopy);
  bool res;

  // Solve J*a = input_x
  res = res & grpPtr->applyRightPreconditioning(params, input_x, *a);

  // Solve J*b = dR/dp
  res = res & grpPtr->applyRightPreconditioning(params, *derivResidualParamPtr, *b);

  // Compute (dJy/dx)*a
  res = res & derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), *a, 
				      tpFVec.getNullVec(), *derivJ);

  // Compute input_null - (dJy/dx)*a
  *derivJ = derivJ->update(1.0, input_null, -1.0);

  // Solve J*c = input_null - (dJy/dx)*a
  res = res & grpPtr->applyRightPreconditioning(params, *derivJ, *c);

  // Compute (dJy/dx)*b
  res = res & derivPtr->computeDJnDxa(*grpPtr, tpXVec.getNullVec(), *b, 
				      tpFVec.getNullVec(), *derivJ);

  // Compute (dJy/dx)*b - dJy/dp
  *derivJ = derivJ->update(-1.0, *derivNullResidualParamPtr, 1.0);

  // Solve J*d = (dJy/dx)*b - dJy/dp
  res = res & grpPtr->applyRightPreconditioning(params, *derivJ, *d);

  // Compute result_param = -(input_param - phi^T*c)/(phi^T*d)
  result_param = (input_param- lengthVecPtr->dot(*c))/(lengthVecPtr->dot(*d));

  // Compute result_x = a - result_param*b 
  *result_x = result_x->update(1.0, *a, -result_param, *b, 0.0);

  // Compute result_null = c + result_param*d
  *result_null = result_null->update(1.0, *c, result_param, *d, 0.0);

  tp_result.setVec(*result_x, *result_null, result_param);

  // Clean up memory
  delete result_x;
  delete result_null;
  delete tmp;

  delete a;
  delete b;
  delete c;
  delete d;
  delete derivJ;

  return res;
}

bool
LOCA::Bifurcation::TPBordGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::Bifurcation::TPBordGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Bifurcation::TPBordGroup::isGradient() const 
{
  return false;
}

bool
LOCA::Bifurcation::TPBordGroup::isNewton() const 
{
  return isValidNewton;
}
  
const Vector&
LOCA::Bifurcation::TPBordGroup::getX() const 
{
  return tpXVec;
}

const Vector&
LOCA::Bifurcation::TPBordGroup::getF() const 
{
  return tpFVec;
}

double
LOCA::Bifurcation::TPBordGroup::getNormF() const 
{
  return tpFVec.norm();
}

const Vector&
LOCA::Bifurcation::TPBordGroup::getGradient() const 
{
  cout << "ERROR: LOCA::Bifurcation::TPBordGroup::getGradient() "
       << " - not implemented" << endl;
  throw "LOCA Error";
  return getNewton();
}

const Vector&
LOCA::Bifurcation::TPBordGroup::getTangent() const 
{
  return tpTangentVec;
}

const Vector&
LOCA::Bifurcation::TPBordGroup::getNewton() const 
{
  return tpNewtonVec;
}

double
LOCA::Bifurcation::TPBordGroup::getNormNewtonSolveResidual() const 
{
  TPBordVector residual = tpFVec;
  
  applyJacobian(tpNewtonVec, residual);
  residual = residual.update(1.0, tpFVec, 1.0);
  return residual.norm();
}

bool
LOCA::Bifurcation::TPBordGroup::print() const 
{
  cout << "Beginning TPBordGroup.print:" << endl;
  cout << "Group = " << endl;
  grpPtr->print();

  cout << endl;

  cout << "tpXVec = " << endl;
  tpXVec.print();

  cout << endl;

  if (isValidF) {
    cout << "tpFVec = " << endl;
    tpFVec.print();
  }
  else
    cout << "tpFVec not computed" << endl;

  cout << endl;

  if (isValidNewton) {
    cout << "tpNewtonVec = " << endl;
    tpNewtonVec.print();
  }
  else
    cout << "tpNewtonVec not computed" << endl;

  cout << endl;

  cout << "lengthVec = " << endl;
  lengthVecPtr->print();

  cout << endl;

  cout << "bifParmId = " << bifParamId << endl;

  cout << endl;

  if (isValidJacobian) {
    cout << "derivResdiualParam = " << endl;
    derivResidualParamPtr->print();
    cout << endl;
    cout << "derivNullResidualParam = " << endl;
    derivNullResidualParamPtr->print();
  }
  else
    cout << "deriv vectors not computed" << endl;

  cout << "End TPBordGroup.print:" << endl;
  cout << endl;

  return true;
}

