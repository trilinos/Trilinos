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

#include "LOCA_Bifurcation_HopfBord_ExtendedGroup.H"
#include "LOCA_Bifurcation_HopfBord_AbstractGroup.H" 
#include "LOCA_Parameter_Vector.H"
#include "NOX_Parameter_List.H"

LOCA::Bifurcation::HopfBord::ExtendedGroup::ExtendedGroup(
			      LOCA::Bifurcation::HopfBord::AbstractGroup& g,
			      const NOX::Abstract::Vector& realEigenVec,
			      const NOX::Abstract::Vector& imaginaryEigenVec,
			      NOX::Abstract::Vector& lenVec,
			      double frequency,
			      int paramId)
  : grpPtr(&g), 
    hopfXVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    hopfFVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    hopfNewtonVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    lengthVecPtr(&lenVec), 
    bifParamId(paramId), 
    derivResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    derivRealEigenResidualParamPtr(lenVec.clone(NOX::ShapeCopy)),
    derivImagEigenResidualParamPtr(lenVec.clone(NOX::ShapeCopy)),
    massTimesYPtr(lenVec.clone(NOX::ShapeCopy)),
    minusMassTimesZPtr(lenVec.clone(NOX::ShapeCopy)),
    ownsGroup(false),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  init();
}

LOCA::Bifurcation::HopfBord::ExtendedGroup::ExtendedGroup(
			   const LOCA::Bifurcation::HopfBord::AbstractGroup& g,
			   const NOX::Abstract::Vector& realEigenVec,
			   const NOX::Abstract::Vector& imaginaryEigenVec,
			   NOX::Abstract::Vector& lenVec,
			   double frequency,
			   int paramId)
  : grpPtr(dynamic_cast<LOCA::Bifurcation::HopfBord::AbstractGroup*>(g.clone())), 
    hopfXVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    hopfFVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    hopfNewtonVec(g.getX(), realEigenVec, imaginaryEigenVec, frequency, 0.0),
    lengthVecPtr(&lenVec), 
    bifParamId(paramId), 
    derivResidualParamPtr(lenVec.clone(NOX::ShapeCopy)), 
    derivRealEigenResidualParamPtr(lenVec.clone(NOX::ShapeCopy)),
    derivImagEigenResidualParamPtr(lenVec.clone(NOX::ShapeCopy)),
    massTimesYPtr(lenVec.clone(NOX::ShapeCopy)),
    minusMassTimesZPtr(lenVec.clone(NOX::ShapeCopy)),
    ownsGroup(true),
    isValidF(false),
    isValidJacobian(false),
    isValidNewton(false)
{
  init();
}

LOCA::Bifurcation::HopfBord::ExtendedGroup::ExtendedGroup(
		    const LOCA::Bifurcation::HopfBord::ExtendedGroup& source, 
		    NOX::CopyType type)
  : grpPtr(dynamic_cast<LOCA::Bifurcation::HopfBord::AbstractGroup*>(source.grpPtr->clone())), 
    hopfXVec(source.hopfXVec, type),
    hopfFVec(source.hopfFVec, type),
    hopfNewtonVec(source.hopfNewtonVec, type),
    lengthVecPtr(source.lengthVecPtr), 
    bifParamId(source.bifParamId),
    derivResidualParamPtr(source.derivResidualParamPtr->clone(type)),
    derivRealEigenResidualParamPtr(source.derivRealEigenResidualParamPtr->clone(type)),
    derivImagEigenResidualParamPtr(source.derivImagEigenResidualParamPtr->clone(type)),
    massTimesYPtr(source.massTimesYPtr->clone(type)),
    minusMassTimesZPtr(source.minusMassTimesZPtr->clone(type)),
    ownsGroup(true),
    isValidF(source.isValidF),
    isValidJacobian(source.isValidJacobian),
    isValidNewton(source.isValidNewton) {}


LOCA::Bifurcation::HopfBord::ExtendedGroup::~ExtendedGroup() 
{
  if (ownsGroup)
    delete grpPtr;
  delete derivResidualParamPtr;
  delete derivRealEigenResidualParamPtr;
  delete derivImagEigenResidualParamPtr;
  delete massTimesYPtr;
  delete minusMassTimesZPtr;
}

LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::HopfBord::ExtendedGroup::operator=(
			     const LOCA::Continuation::AbstractGroup& source)
{
  return *this = dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedGroup&>(source);
}

NOX::Abstract::Group&
LOCA::Bifurcation::HopfBord::ExtendedGroup::operator=(
				         const NOX::Abstract::Group& source)
{
  return *this = dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedGroup&>(source);
}

LOCA::Bifurcation::HopfBord::ExtendedGroup&
LOCA::Bifurcation::HopfBord::ExtendedGroup::operator=(
			  const LOCA::Bifurcation::HopfBord::ExtendedGroup& source) 
{

  // Protect against A = A
  if (this != &source) {
    NOX::CopyType type = NOX::DeepCopy;

    // Delete old values
    delete derivResidualParamPtr;
    delete derivRealEigenResidualParamPtr;
    delete derivImagEigenResidualParamPtr;

    // Copy values
    *grpPtr = *source.grpPtr;
    hopfXVec = source.hopfXVec;
    hopfFVec = source.hopfFVec;
    hopfNewtonVec = source.hopfNewtonVec;
    lengthVecPtr = source.lengthVecPtr;
    derivResidualParamPtr = source.derivResidualParamPtr->clone(type);
    derivRealEigenResidualParamPtr = source.derivRealEigenResidualParamPtr->clone(type);
    derivImagEigenResidualParamPtr = source.derivImagEigenResidualParamPtr->clone(type);
    massTimesYPtr = source.massTimesYPtr->clone(type);
    minusMassTimesZPtr = source.minusMassTimesZPtr->clone(type);
    bifParamId = source.bifParamId;
    isValidF = source.isValidF;
    isValidJacobian = source.isValidJacobian;
    isValidNewton = source.isValidNewton;
  }

  return *this;
}

NOX::Abstract::Group*
LOCA::Bifurcation::HopfBord::ExtendedGroup::clone(NOX::CopyType type) const 
{
  return new LOCA::Bifurcation::HopfBord::ExtendedGroup(*this, type);
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setParams(const LOCA::ParameterVector& p) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(p);
}

const LOCA::ParameterVector&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getParams() const 
{
  return grpPtr->getParams();
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setParam(int paramID, double val)
{
  grpPtr->setParam(paramID, val);
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::getParam(int paramID) const
{
  return grpPtr->getParam(paramID);
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setParam(string paramID, double val)
{
  grpPtr->setParam(paramID, val);
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::getParam(string paramID) const
{
  return grpPtr->getParam(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeDfDp(int paramID, 
					      NOX::Abstract::Vector& result)
{
  NOX::Abstract::Group::ReturnType res;

  // Cast result to Hopf vector
  LOCA::Bifurcation::HopfBord::ExtendedVector& hopf_result = 
    dynamic_cast<LOCA::Bifurcation::HopfBord::ExtendedVector&>(result);

  // Compute f, J, (J+iwB)*(y+iz)
  res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute df/dp
  res = grpPtr->computeDfDp(paramID, hopf_result.getXVec());
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute [(J+iwB)*(y+iz)]_p
  res = grpPtr->computeDCeDp(hopfXVec.getRealEigenVec(),
			     hopfXVec.getImagEigenVec(),
			     hopfXVec.getFrequency(),
			     paramID,
			     hopfFVec.getRealEigenVec(),
			     hopfFVec.getImagEigenVec(),
			     hopf_result.getRealEigenVec(),
			     hopf_result.getImagEigenVec());
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Set frequency  component
  hopf_result.getFrequency() = 0.0;

  // Set parameter component
  hopf_result.getBifParam() = 0.0;

  return res;
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::getBifParam() const 
{
  LOCA::ParameterVector params(grpPtr->getParams());
  return params[bifParamId];
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setBifParam(double param) 
{
  LOCA::ParameterVector params(grpPtr->getParams());

  params[bifParamId] = param;
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->setParams(params);
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setX(const NOX::Abstract::Vector& y) 
{
  setX( dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector&>(y) );
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::setX(
			       const LOCA::Bifurcation::HopfBord::ExtendedVector& y) 
{
  grpPtr->setX( y.getXVec() );
  hopfXVec = y;

  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeX(
					     const NOX::Abstract::Group& g, 
					     const NOX::Abstract::Vector& d,
					     double step) 
{
  computeX( dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedGroup&>(g),
	    dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector&>(d),
	    step);
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeX(
			 const LOCA::Bifurcation::HopfBord::ExtendedGroup& g, 
			 const LOCA::Bifurcation::HopfBord::ExtendedVector& d,
			 double step) 
{
  isValidF = false;
  isValidJacobian = false;
  isValidNewton = false;

  grpPtr->computeX(*(g.grpPtr), d.getXVec(), step);
  hopfXVec.update(1.0, g.getX(), step, d, 0.0);
  setBifParam(hopfXVec.getBifParam());
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeF() 
{
  if (isValidF)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res;

  // Compute F
  res = grpPtr->computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  hopfFVec.getXVec() = grpPtr->getF();
  
  // Compute J
  res = grpPtr->computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute B
  res = grpPtr->computeMassMatrix();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute (J+iwB)*(y+iz)
  res = grpPtr->applyComplex(hopfXVec.getRealEigenVec(),
			     hopfXVec.getImagEigenVec(),
			     hopfXVec.getFrequency(),
			     hopfFVec.getRealEigenVec(),
			     hopfFVec.getImagEigenVec());
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  // Compute l^T*y - 1
  hopfFVec.getFrequency() = hopfXVec.getRealEigenVec().dot(*lengthVecPtr)-1.0;

  // Compute l^T*z
  hopfFVec.getBifParam() = hopfXVec.getImagEigenVec().dot(*lengthVecPtr);
  
  isValidF = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeJacobian() 
{
  if (isValidJacobian)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute J
  res = grpPtr->computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute B
  res = grpPtr->computeMassMatrix();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  // Compute dF/dp
  res = grpPtr->computeDfDp(bifParamId, *derivResidualParamPtr);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute [(J+iwB)*(y+iz)]_p
  res = grpPtr->computeDCeDp(hopfXVec.getRealEigenVec(),
			     hopfXVec.getImagEigenVec(),
			     hopfXVec.getFrequency(),
			     bifParamId,
			     hopfFVec.getRealEigenVec(),
			     hopfFVec.getImagEigenVec(),
			     *derivRealEigenResidualParamPtr,
			     *derivImagEigenResidualParamPtr);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute B*y
  res = grpPtr->applyMassMatrix(hopfXVec.getRealEigenVec(),
				*massTimesYPtr);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute -B*z
  res = grpPtr->applyMassMatrix(hopfXVec.getImagEigenVec(),
				*minusMassTimesZPtr);
  if (res != NOX::Abstract::Group::Ok)
    return res;
  minusMassTimesZPtr->scale(-1.0);

  isValidJacobian = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeGradient() 
{
  return NOX::Abstract::Group::NotDefined;
}
   
NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::computeNewton(
						NOX::Parameter::List& params) 
{
  if (isValidNewton)
    return NOX::Abstract::Group::Ok;

  NOX::Abstract::Group::ReturnType res = computeF();
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  res = computeJacobian();
  if (res != NOX::Abstract::Group::Ok)
    return res;

  res = applyJacobianInverse(params, hopfFVec, hopfNewtonVec);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  hopfNewtonVec.scale(-1.0);
  isValidNewton = true;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::applyJacobian(
					  const NOX::Abstract::Vector& input,
					  NOX::Abstract::Vector& result) const 
{
  // Cast vectors to HopfBordVectors
  const LOCA::Bifurcation::HopfBord::ExtendedVector& hopf_input = 
    dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector&>(input);
  LOCA::Bifurcation::HopfBord::ExtendedVector& hopf_result = 
    dynamic_cast<LOCA::Bifurcation::HopfBord::ExtendedVector&>(result);

  // Get constant references to input vector components
  const NOX::Abstract::Vector& input_x = hopf_input.getXVec();
  const NOX::Abstract::Vector& input_y = hopf_input.getRealEigenVec();
  const NOX::Abstract::Vector& input_z = hopf_input.getImagEigenVec();
  double input_w = hopf_input.getFrequency();
  double input_param = hopf_input.getBifParam();

  // Get non-constant references to result vector components
  NOX::Abstract::Vector& result_x = hopf_result.getXVec();
  NOX::Abstract::Vector& result_y = hopf_result.getRealEigenVec();
  NOX::Abstract::Vector& result_z = hopf_result.getImagEigenVec();
  double& result_w = hopf_result.getFrequency();
  double& result_param = hopf_result.getBifParam();

  // Temporary vectors
  NOX::Abstract::Vector *tmp_real = input_x.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *tmp_imag = input_x.clone(NOX::ShapeCopy);

  // Return value
  NOX::Abstract::Group::ReturnType res;

  // compute J*X
  res = grpPtr->applyJacobian(input_x, result_x);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // compute J*X + P*dR/dp
  result_x.update(input_param, *derivResidualParamPtr, 1.0);

  // compute (J+iwB)*(Y+iZ)
  res = grpPtr->applyComplex(input_y, input_z, hopfXVec.getFrequency(), 
			     result_y, result_z);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // compute (J+iwB)*(Y+iZ) + P*d(J+iwB)*(y+iz)/dp + iW*B*(y+iz)
  result_y.update(input_param, *derivRealEigenResidualParamPtr, 
		  input_w, *minusMassTimesZPtr, 1.0);
  result_z.update(input_param, *derivImagEigenResidualParamPtr, 
		  input_w, *massTimesYPtr, 1.0);

  // compute (d(J+iwB)*(y+iz)/dx)*X
  res = grpPtr->computeDCeDxa(hopfXVec.getRealEigenVec(), 
			      hopfXVec.getImagEigenVec(), 
			      hopfXVec.getFrequency(),
			      input_x, 
			      hopfFVec.getRealEigenVec(), 
			      hopfFVec.getImagEigenVec(),
			      *tmp_real,
			      *tmp_imag);
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  // (d(J+iwB)*(y+iz)/dx)*X + (J+iwB)*(Y+iZ) + P*d(J+iwB)*(y+iz)/dp + iW*B*(y+iz)
  result_y.update(1.0, *tmp_real, 1.0);
  result_z.update(1.0, *tmp_imag, 1.0);

  // compute phi^T*Y, phi^T*Z
  result_w = lengthVecPtr->dot(input_y);
  result_param = lengthVecPtr->dot(input_z);

  delete tmp_real;
  delete tmp_imag;

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::applyJacobianTranspose(
					 const NOX::Abstract::Vector& input,
					 NOX::Abstract::Vector& result) const 
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::applyJacobianInverse(
					 NOX::Parameter::List& params,
					 const NOX::Abstract::Vector& input,
					 NOX::Abstract::Vector& result) const 
{
  const NOX::Abstract::Vector* inputs[1];
  NOX::Abstract::Vector* results[1];

  inputs[0] = &input;
  results[0] = &result;

  return applyJacobianInverseMulti(params, inputs, results, 1);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::ExtendedGroup::applyJacobianInverseMulti(
			    NOX::Parameter::List& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** results, int nVecs) const 
{
  // Number of input vectors
  int m = nVecs; 

  // Return type
   NOX::Abstract::Group::ReturnType res;
  
  // Build arrays of solution, eigenvector and parameter components
  const NOX::Abstract::Vector** inputs_x = 
    new const NOX::Abstract::Vector*[m+1];
  const NOX::Abstract::Vector** inputs_y =
    new const NOX::Abstract::Vector*[m+1];
  const NOX::Abstract::Vector** inputs_z =
    new const NOX::Abstract::Vector*[m+1];
  double* inputs_w = new double[m];
  double* inputs_params = new double[m];

  NOX::Abstract::Vector** results_x = new NOX::Abstract::Vector*[m+1];
  NOX::Abstract::Vector** results_y = new NOX::Abstract::Vector*[m+2];
  NOX::Abstract::Vector** results_z = new NOX::Abstract::Vector*[m+2];
  double** results_w = new double*[m];
  double** results_params = new double*[m];

  NOX::Abstract::Vector** tmp_real = new NOX::Abstract::Vector*[m+2];
  NOX::Abstract::Vector** tmp_imag = new NOX::Abstract::Vector*[m+2];

  const LOCA::Bifurcation::HopfBord::ExtendedVector* constHopfVecPtr;
  LOCA::Bifurcation::HopfBord::ExtendedVector* hopfVecPtr;

  for (int i=0; i<m; i++) {
    constHopfVecPtr = 
      dynamic_cast<const LOCA::Bifurcation::HopfBord::ExtendedVector*>(inputs[i]);
    inputs_x[i] = &(constHopfVecPtr->getXVec());
    inputs_y[i] = &(constHopfVecPtr->getRealEigenVec());
    inputs_z[i] = &(constHopfVecPtr->getImagEigenVec());
    inputs_w[i] = constHopfVecPtr->getFrequency();
    inputs_params[i] = constHopfVecPtr->getBifParam();
 
    hopfVecPtr = 
      dynamic_cast<LOCA::Bifurcation::HopfBord::ExtendedVector*>(results[i]);
    results_x[i] = &(hopfVecPtr->getXVec());
    results_y[i] = &(hopfVecPtr->getRealEigenVec());
    results_z[i] = &(hopfVecPtr->getImagEigenVec());
    results_w[i] = &(hopfVecPtr->getFrequency());
    results_params[i] = &(hopfVecPtr->getBifParam());

    tmp_real[i] = inputs_x[i]->clone(NOX::ShapeCopy);
    tmp_imag[i] = inputs_x[i]->clone(NOX::ShapeCopy);
  }

  // Set next to last components to deriv. w.r.t. parameter
  inputs_x[m] = derivResidualParamPtr;
  inputs_y[m] = derivRealEigenResidualParamPtr;
  inputs_z[m] = derivImagEigenResidualParamPtr;
  results_x[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  results_y[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  results_z[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  tmp_real[m] = inputs_x[m]->clone(NOX::ShapeCopy);
  tmp_imag[m] = inputs_x[m]->clone(NOX::ShapeCopy);

  results_y[m+1] = inputs_y[m+1]->clone(NOX::ShapeCopy);
  results_z[m+1] = inputs_y[m+1]->clone(NOX::ShapeCopy);
  tmp_real[m+1] = inputs_x[m+1]->clone(NOX::ShapeCopy);
  tmp_imag[m+1] = inputs_x[m+1]->clone(NOX::ShapeCopy);

  // Solve J*results_x = inputs_x
  res = applyJacobianInverseMulti(params, inputs_x, results_x, m+1);
  if (res != NOX::Abstract::Group::Ok)
    return res;
  

  // Compute tmp = (inputs_y + i*inputs_z) - (dJ(J+iwB)(y+iz)/dx)*results_x
  for (int i=0; i<m+1; i++) {
    res = grpPtr->computeDCeDxa(hopfXVec.getRealEigenVec(), 
				hopfXVec.getImagEigenVec(),
				hopfXVec.getFrequency(),
				*results_x[i],
				hopfFVec.getRealEigenVec(), 
				hopfFVec.getImagEigenVec(),
				*tmp_real[i],
				*tmp_imag[i]);
    if (res != NOX::Abstract::Group::Ok)
      return res;

    tmp_real[i]->update(1.0, *inputs_y[i], -1.0);
    tmp_imag[i]->update(1.0, *inputs_z[i], -1.0);
  } 
  tmp_real[m+1] = minusMassTimesZPtr;
  tmp_imag[m+1] = massTimesYPtr;

  // Solve (J+iwB)*(results_y + i*results_z) = (tmp_real + i*tmp_imag)
  res = grpPtr->applyComplexInverseMulti(params, tmp_real, tmp_imag, 
					 hopfXVec.getFrequency(),
					 results_y, results_z, m+2);
  if (res != NOX::Abstract::Group::Ok)
    return res;
  
  // Compute and set results
  double ltc = lengthVecPtr->dot(*results_y[m+1]);
  double ltd = lengthVecPtr->dot(*results_z[m+1]);
  double ltg = lengthVecPtr->dot(*results_y[m]);
  double lth = lengthVecPtr->dot(*results_z[m]);
  double denom_param = ltd*ltg - ltc*lth;
  double lte, ltf;
 
  for (int i=0; i<m; i++) {
    lte = lengthVecPtr->dot(*results_y[i]);
    ltf = lengthVecPtr->dot(*results_z[i]);
    *results_params[i] = 
      (ltc*(inputs_params[i] - ltf) + ltd*(lte - inputs_w[i])) / denom_param;
    *results_w[i] = (inputs_params[i] - ltf + lth*(*results_params[i]))/ltd;

    results_x[i]->update(-*results_params[i], *results_x[m], 1.0);
    results_y[i]->update(-*results_params[i], *results_y[m], 
			 -*results_w[i], *results_y[m+1], 1.0);
    results_z[i]->update(-*results_params[i], *results_z[m], 
			 -*results_w[i], *results_z[m+1], 1.0);

    delete tmp_real[i];
    delete tmp_imag[i];
  }

  delete results_x[m];
  delete results_y[m];
  delete results_z[m];
  delete tmp_real[m];
  delete tmp_imag[m];

  delete results_y[m+1];
  delete results_z[m+1];

  delete [] inputs_x;
  delete [] inputs_y;
  delete [] inputs_z;
  delete [] inputs_w;
  delete [] inputs_params;
  delete [] results_x;
  delete [] results_y;
  delete [] results_z;
  delete [] results_w;
  delete [] results_params;
  delete [] tmp_real;
  delete [] tmp_imag;

  return res;
}

bool
LOCA::Bifurcation::HopfBord::ExtendedGroup::isF() const 
{
  return isValidF;
}

bool
LOCA::Bifurcation::HopfBord::ExtendedGroup::isJacobian() const 
{
  return isValidJacobian;
}

bool
LOCA::Bifurcation::HopfBord::ExtendedGroup::isGradient() const 
{
  return false;
}

bool
LOCA::Bifurcation::HopfBord::ExtendedGroup::isNewton() const 
{
  return isValidNewton;
}
  
const NOX::Abstract::Vector&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getX() const 
{
  return hopfXVec;
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getF() const 
{
  return hopfFVec;
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::getNormF() const 
{
  return hopfFVec.norm();
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getGradient() const 
{
  cout << "ERROR: LOCA::Bifurcation::HopfBord::ExtendedGroup::getGradient() "
       << " - not implemented" << endl;
  throw "LOCA Error";
  return getNewton();
}

const NOX::Abstract::Vector&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getNewton() const 
{
  return hopfNewtonVec;
}

double
LOCA::Bifurcation::HopfBord::ExtendedGroup::getNormNewtonSolveResidual() const 
{
  LOCA::Bifurcation::HopfBord::ExtendedVector residual = hopfFVec;
  
  NOX::Abstract::Group::ReturnType res = applyJacobian(hopfNewtonVec,residual);
  if (res != NOX::Abstract::Group::Ok) {
    cout << "ERROR: applyJacobian() in getNormNewtonSolveResidual "
	 << " returned not ok" << endl;
    throw "LOCA Error";
    return 0.0;
  }

  residual.update(1.0, hopfFVec, 1.0);
  return residual.norm();
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::printSolution(const double conParam) const 
{
  cout << "LOCA::Bifurcation::HopfBord::ExtendedGroup::printSolution\n";

  cout << "\tPrinting Solution Vector for conParam = " << conParam << endl;
  grpPtr->printSolution(conParam);

  cout << "\tPrinting Real Component of Eigenvector for bif param = " 
       << getBifParam() << endl;
  grpPtr->printSolution(hopfXVec.getRealEigenVec(), hopfXVec.getBifParam());

  cout << "\tPrinting Imaginary Component of Eigenvector for bif param = " 
       << getBifParam() << endl;
  grpPtr->printSolution(hopfXVec.getImagEigenVec(), 
			hopfXVec.getBifParam());
}

const LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getUnderlyingGroup() const
{
  return grpPtr->getUnderlyingGroup();
}

LOCA::Continuation::AbstractGroup&
LOCA::Bifurcation::HopfBord::ExtendedGroup::getUnderlyingGroup()
{
  return grpPtr->getUnderlyingGroup();
}

void
LOCA::Bifurcation::HopfBord::ExtendedGroup::init()
{
  hopfXVec.getBifParam() = getBifParam();

  // Rescale length vector so that the normalization condition is met
  double lVecDotY;
  lVecDotY = hopfXVec.getRealEigenVec().dot(*lengthVecPtr);
  if (lVecDotY == 0.0) {
    cout << "ERROR: LOCA::Bifurcation::HopfBord::ExtendedGroup::HopfPBordGroup\n"
         << "     : length vector can not have Norm zero " << endl;

    throw "LOCA Error";
  }
  cout << "\tIn HopfBord::Group Constructor, scaling lenVec by:" 
       << (1.0 / lVecDotY) << endl;
  lengthVecPtr->scale(1.0 / lVecDotY);
}
