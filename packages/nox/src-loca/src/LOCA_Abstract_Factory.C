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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "LOCA_Abstract_Factory.H"

bool
LOCA::Abstract::Factory::createPredictorStrategy(
        const string& strategyName,
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& predictorParams,
	Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createContinuationStrategy(
    const string& strategyName,
    const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RefCountPtr<Teuchos::ParameterList>& stepperParams,
    const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
    const Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>& pred,
    const vector<int>& paramIDs,
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createBifurcationStrategy(
    const string& strategyName,
    const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RefCountPtr<Teuchos::ParameterList>& bifurcationParams,
    const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createStepSizeStrategy(
        const string& strategyName,
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& stepsizeParams,
	Teuchos::RefCountPtr<LOCA::StepSize::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createBorderedSolverStrategy(
        const string& strategyName,
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams,
	Teuchos::RefCountPtr<LOCA::BorderedSolver::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createEigensolverStrategy(
         const string& strategyName,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams,
	 Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createEigenvalueSortStrategy(
        const string& strategyName,
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams,
	Teuchos::RefCountPtr<LOCA::EigenvalueSort::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createSaveEigenDataStrategy(
         const string& strategyName,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams,
	 Teuchos::RefCountPtr<LOCA::SaveEigenData::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createAnasaziOperatorStrategy(
      const string& strategyName,
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams,
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams,
      const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp,
      Teuchos::RefCountPtr<LOCA::AnasaziOperator::AbstractStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createMooreSpenceTurningPointSolverStrategy(
       const string& strategyName,
       const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams,
       Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::SolverStrategy>& strategy)
{
  return false;
}

bool
LOCA::Abstract::Factory::createMooreSpencePitchforkSolverStrategy(
       const string& strategyName,
       const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams,
       Teuchos::RefCountPtr<LOCA::Pitchfork::MooreSpence::SolverStrategy>& strategy)
{
  return false;
}
