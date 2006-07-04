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

#include "LOCA_MultiPredictor_Constant.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_MultiContinuation_ExtendedMultiVector.H"

LOCA::MultiPredictor::Constant::Constant(
	      const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	      const Teuchos::RefCountPtr<Teuchos::ParameterList>& predParams) :
  globalData(global_data),
  predictor(),
  secant(),
  initialized(false)
{
}

LOCA::MultiPredictor::Constant::~Constant()
{
}

LOCA::MultiPredictor::Constant::Constant(
				 const LOCA::MultiPredictor::Constant& source,
				 NOX::CopyType type) :
  globalData(source.globalData),
  predictor(),
  secant(),
  initialized(source.initialized)
{
  if (source.initialized) {
    predictor = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(source.predictor->clone(type));

    secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(source.secant->clone(type));
  }
}

LOCA::MultiPredictor::AbstractStrategy&
LOCA::MultiPredictor::Constant::operator=(
			  const LOCA::MultiPredictor::AbstractStrategy& s)
{
  const LOCA::MultiPredictor::Constant& source = 
    dynamic_cast<const LOCA::MultiPredictor::Constant&>(s);

  if (this != &source) {
    globalData = source.globalData;
    initialized = source.initialized;

    if (source.initialized) {
      predictor = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(source.predictor->clone(NOX::DeepCopy));

      secant = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(source.secant->clone(NOX::DeepCopy));
    }
  }

  return *this;
}

Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>
LOCA::MultiPredictor::Constant::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Constant(*this, type));
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Constant::compute(
	      bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      const LOCA::MultiContinuation::ExtendedVector& prevXVec,
	      const LOCA::MultiContinuation::ExtendedVector& xVec)
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails))
    globalData->locaUtils->out() << 
      "\n\tCalling Predictor with method: Constant" << std::endl;

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
  for (int i=0; i<numParams; i++)
    predictor->getScalar(i,i) = 1.0;

  // Set orientation based on parameter change
  setPredictorOrientation(baseOnSecant, stepSize, grp, prevXVec, 
			  xVec, *secant, *predictor);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Constant::evaluate(
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
LOCA::MultiPredictor::Constant::computeTangent(
			LOCA::MultiContinuation::ExtendedMultiVector& v)
{
  v = *predictor;

  return NOX::Abstract::Group::Ok;
}

bool
LOCA::MultiPredictor::Constant::isTangentScalable() const
{
  return false;
}
