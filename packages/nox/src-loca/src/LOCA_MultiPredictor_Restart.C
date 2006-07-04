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
#include "LOCA_MultiPredictor_Restart.H"
#include "LOCA_GlobalData.H"
#include "NOX_Utils.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_MultiContinuation_ExtendedMultiVector.H"

LOCA::MultiPredictor::Restart::Restart(
	      const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	      const Teuchos::RefCountPtr<Teuchos::ParameterList>& predParams) :
  globalData(global_data),
  predictor()
{
  const char *func = "LOCA::MultiPredictor::Restart::Restart()";

  // Get predictor vector from parameter list
  string name = "Restart Vector";
  if (!predParams->isParameter(name))
    globalData->locaErrorCheck->throwError(func, name + " is not set!");

  if ((*predParams).INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedMultiVector> >(name)) 
    predictor = (*predParams).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedMultiVector> >(name);

  else if ((*predParams).INVALID_TEMPLATE_QUALIFIER
	   isType< Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedVector> >(name)) {
    Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedVector> v =
      (*predParams).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedVector> >(name);
    predictor = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(v->createMultiVector(1, NOX::DeepCopy));
  }
  else
    globalData->locaErrorCheck->throwError(func, name + " is not a Teuchos::RefCountPtr to a LOCA::Extended::Vector nor a LOCA::Extended::MultiVector!");

  // Note we don't need a secant vector since it is assumed the orientation
  // is already correct
}

LOCA::MultiPredictor::Restart::~Restart()
{
}

LOCA::MultiPredictor::Restart::Restart(
				 const LOCA::MultiPredictor::Restart& source,
				 NOX::CopyType type) :
  globalData(source.globalData),
  predictor(source.predictor)
{
}

LOCA::MultiPredictor::AbstractStrategy&
LOCA::MultiPredictor::Restart::operator=(
			  const LOCA::MultiPredictor::AbstractStrategy& s)
{
  const LOCA::MultiPredictor::Restart& source = 
    dynamic_cast<const LOCA::MultiPredictor::Restart&>(s);

  if (this != &source) {
    globalData = source.globalData;
    predictor = source.predictor;
  }

  return *this;
}

Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>
LOCA::MultiPredictor::Restart::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Restart(*this, type));
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Restart::compute(
	      bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      const LOCA::MultiContinuation::ExtendedVector& prevXVec,
	      const LOCA::MultiContinuation::ExtendedVector& xVec)
{
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperDetails))
    globalData->locaUtils->out() << 
      "\n\tCalling Predictor with method: Restart" << std::endl;

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Restart::evaluate(
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
LOCA::MultiPredictor::Restart::computeTangent(
			LOCA::MultiContinuation::ExtendedMultiVector& v)
{
  v = *predictor;

  return NOX::Abstract::Group::Ok;
}

bool
LOCA::MultiPredictor::Restart::isTangentScalable() const
{
  return false;
}
