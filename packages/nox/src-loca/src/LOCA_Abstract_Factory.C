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

#include "LOCA_Abstract_Factory.H"

bool
LOCA::Abstract::Factory::createPredictorStrategy(
        const std::string& strategyName,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& predictorParams,
	Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createContinuationStrategy(
    const std::string& strategyName,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& stepperParams,
    const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
    const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
    const std::vector<int>& paramIDs,
    Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createBifurcationStrategy(
    const std::string& strategyName,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& bifurcationParams,
    const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createStepSizeStrategy(
        const std::string& strategyName,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams,
	Teuchos::RCP<LOCA::StepSize::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createBorderedSolverStrategy(
        const std::string& strategyName,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
	Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createEigensolverStrategy(
         const std::string& strategyName,
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
	 Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createEigenvalueSortStrategy(
        const std::string& strategyName,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
	Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createSaveEigenDataStrategy(
         const std::string& strategyName,
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
	 Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createAnasaziOperatorStrategy(
      const std::string& strategyName,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
      const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
      const Teuchos::RCP<NOX::Abstract::Group>& grp,
      Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createMooreSpenceTurningPointSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createMooreSpencePitchforkSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::Pitchfork::MooreSpence::SolverStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createMooreSpenceHopfSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::Hopf::MooreSpence::SolverStrategy>& strategy)
{
  return false;
}
