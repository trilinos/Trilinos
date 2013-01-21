// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
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
	      const Teuchos::RCP<LOCA::GlobalData>& global_data,
	      const Teuchos::RCP<Teuchos::ParameterList>& predParams) :
  globalData(global_data),
  predictor()
{
  const char *func = "LOCA::MultiPredictor::Restart::Restart()";

  // Get predictor vector from parameter list
  std::string name = "Restart Vector";
  if (!predParams->isParameter(name))
    globalData->locaErrorCheck->throwError(func, name + " is not set!");

  if ((*predParams).INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RCP<LOCA::MultiContinuation::ExtendedMultiVector> >(name)) 
    predictor = (*predParams).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<LOCA::MultiContinuation::ExtendedMultiVector> >(name);

  else if ((*predParams).INVALID_TEMPLATE_QUALIFIER
	   isType< Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> >(name)) {
    Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> v =
      (*predParams).INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> >(name);
    predictor = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedMultiVector>(v->createMultiVector(1, NOX::DeepCopy));
  }
  else
    globalData->locaErrorCheck->throwError(func, name + " is not a Teuchos::RCP to a LOCA::Extended::Vector nor a LOCA::Extended::MultiVector!");

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

Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
LOCA::MultiPredictor::Restart::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Restart(*this, type));
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiPredictor::Restart::compute(
	      bool baseOnSecant, const std::vector<double>& stepSize,
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
	      const std::vector<double>& stepSize,
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
