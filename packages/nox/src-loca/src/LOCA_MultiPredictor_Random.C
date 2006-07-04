// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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

#include "Teuchos_ParameterList.hpp"
#include "LOCA_MultiPredictor_Random.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_MultiContinuation_ExtendedMultiVector.H"

LOCA::MultiPredictor::Random::Random(
	      const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	      const Teuchos::RefCountPtr<Teuchos::ParameterList>& predParams) :
  globalData(global_data),
  predictor(),
  secant(),
  initialized(false),
  epsilon(predParams->get("Epsilon", 1.0e-3))
{
}

LOCA::MultiPredictor::Random::~Random()
{
}

LOCA::MultiPredictor::Random::Random(
				 const LOCA::MultiPredictor::Random& source,
				 NOX::CopyType type) :
  globalData(source.globalData),
  predictor(),
  secant(),
  initialized(source.initialized),
  epsilon(source.epsilon)
{
  if (source.initialized) {
    predictor = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(source.predictor->clone(type));

    secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(source.secant->clone(type));
  }
}

LOCA::MultiPredictor::AbstractStrategy&
LOCA::MultiPredictor::Random::operator=(
			  const LOCA::MultiPredictor::AbstractStrategy& s)
{
  const LOCA::MultiPredictor::Random& source = 
    dynamic_cast<const LOCA::MultiPredictor::Random&>(s);

  if (this != &source) {
    globalData = source.globalData;
    initialized = source.initialized;
    epsilon = source.epsilon;

    if (source.initialized) {
      predictor = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(source.predictor->clone(NOX::DeepCopy));

      secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(source.secant->clone(NOX::DeepCopy));
    }
  }

  return *this;
}

Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>
LOCA::MultiPredictor::Random::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Random(*this, type));
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Random::compute(
	      bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      const LOCA::MultiContinuation::ExtendedVector& prevXVec,
	      const LOCA::MultiContinuation::ExtendedVector& xVec)
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails))
    globalData->locaUtils->out() << 
      "\n\tCalling Predictor with method: Random" << std::endl;

  // Number of continuation parameters
  int numParams = stepSize.size();

  if (!initialized) {

    // Allocate predictor vector
    predictor = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(xVec.createMultiVector(numParams, NOX::ShapeCopy));
    
    // Allocate secant
    secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(xVec.clone(NOX::ShapeCopy));

    initialized = true;
  }

  predictor->init(0.0);

  // Get references to solution components
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> v_x = 
    predictor->getXMultiVec();
  const LOCA::MultiContinuation::ExtendedVector mx = 
    dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(xVec);
  Teuchos::RefCountPtr<const NOX::Abstract::Vector> x_x = mx.getXVec();

  // Fill predictor with random values
  v_x->random();

  for (int i=0; i<numParams; i++) {

    // Scale predictor by solution vector
    (*v_x)[i].scale(*x_x);

    // Scale predictor by epsilon
    (*v_x)[i].scale(epsilon);

  }

  for (int i=0; i<numParams; i++)
    predictor->getScalar(i,i) = 1.0;

  // Set orientation based on parameter change
  setPredictorOrientation(baseOnSecant, stepSize, grp, prevXVec, 
			  xVec, *secant, *predictor);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Random::evaluate(
	      const vector<double>& stepSize,
	      const LOCA::MultiContinuation::ExtendedVector& xVec,
	      LOCA::MultiContinuation::ExtendedMultiVector& result) const
{
  // Number of continuation parameters
  int numParams = stepSize.size();

  for (int i=0; i<numParams; i++)
    result[i].update(1.0, xVec, stepSize[i], (*predictor)[i], 0.0);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Random::computeTangent(
			LOCA::MultiContinuation::ExtendedMultiVector& v)
{
  v = *predictor;

  return NOX::Abstract::Group::Ok;
}

bool
LOCA::MultiPredictor::Random::isTangentScalable() const
{
  return false;
}
