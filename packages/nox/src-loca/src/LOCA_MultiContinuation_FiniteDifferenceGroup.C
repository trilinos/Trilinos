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

#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"

LOCA::MultiContinuation::FiniteDifferenceGroup::FiniteDifferenceGroup()
{
}

LOCA::MultiContinuation::FiniteDifferenceGroup::FiniteDifferenceGroup(
                const LOCA::MultiContinuation::FiniteDifferenceGroup& source, 
		NOX::CopyType type)
{
  if (source.derivPtr != Teuchos::null)
    derivPtr = source.derivPtr->clone(type);
}

LOCA::MultiContinuation::FiniteDifferenceGroup::~FiniteDifferenceGroup() 
{
}

void
LOCA::MultiContinuation::FiniteDifferenceGroup::copy(
					    const NOX::Abstract::Group& src)
{
  const LOCA::MultiContinuation::FiniteDifferenceGroup& source = 
    dynamic_cast<const LOCA::MultiContinuation::FiniteDifferenceGroup&>(src);

  if (this != &source)
    if (source.derivPtr != Teuchos::null)
      derivPtr = source.derivPtr->clone();
}

NOX::Abstract::Group&
LOCA::MultiContinuation::FiniteDifferenceGroup::operator=(
					    const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

void
LOCA::MultiContinuation::FiniteDifferenceGroup::setDerivUtils(
			  const Teuchos::RCP<LOCA::DerivUtils>& deriv)
{
  derivPtr = deriv;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::FiniteDifferenceGroup::computeDfDpMulti(
					  const std::vector<int>& paramIDs, 
					  NOX::Abstract::MultiVector& dfdp, 
					  bool isValidF)
{
  return derivPtr->computeDfDp(*this, paramIDs, dfdp, isValidF);
}

