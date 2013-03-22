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
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_NaturalConstraint.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"

LOCA::MultiContinuation::NaturalGroup::NaturalGroup(
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& continuationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs)
  : LOCA::MultiContinuation::ExtendedGroup(global_data, topParams,
					   continuationParams,
					   grp, pred, paramIDs)
{
  bool skip_dfdp = continuationParams->get("Skip Parameter Derivative", true);
  Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> cons 
    = Teuchos::rcp(new LOCA::MultiContinuation::NaturalConstraint(
	globalData, Teuchos::rcp(this, false)));
  LOCA::MultiContinuation::ExtendedGroup::setConstraints(cons, skip_dfdp);
}

LOCA::MultiContinuation::NaturalGroup::NaturalGroup(
			 const LOCA::MultiContinuation::NaturalGroup& source,
			 NOX::CopyType type)
  : LOCA::MultiContinuation::ExtendedGroup(source, type)
{
  Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::NaturalConstraint>(conGroup->getConstraints())->setNaturalGroup(Teuchos::rcp(this, false));
}


LOCA::MultiContinuation::NaturalGroup::~NaturalGroup() 
{
}

NOX::Abstract::Group&
LOCA::MultiContinuation::NaturalGroup::operator=(
					  const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::MultiContinuation::NaturalGroup::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new LOCA::MultiContinuation::NaturalGroup(*this, type));
}

void
LOCA::MultiContinuation::NaturalGroup::copy(const NOX::Abstract::Group& src) 
{

  // Protect against A = A
  if (this != &src) {
    LOCA::MultiContinuation::ExtendedGroup::copy(src);
  }
}


