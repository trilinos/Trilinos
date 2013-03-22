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

#include "LOCA_Hopf_MooreSpence_FiniteDifferenceGroup.H"

LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup()
{
}

LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::FiniteDifferenceGroup(
            const LOCA::Hopf::MooreSpence::FiniteDifferenceGroup& source, 
	    NOX::CopyType type)
  : LOCA::MultiContinuation::FiniteDifferenceGroup(source, type),
    LOCA::TurningPoint::MooreSpence::FiniteDifferenceGroup(source, type)
{
}


LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::~FiniteDifferenceGroup() 
{
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDp(
			        const std::vector<int>& paramIDs,	     
				const NOX::Abstract::Vector& yVector,
				const NOX::Abstract::Vector& zVector,
				double w,
				NOX::Abstract::MultiVector& result_real,
				NOX::Abstract::MultiVector& result_imag,
				bool isValid)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDCeDp(*this, paramIDs, yVector, zVector, w, 
		 result_real, result_imag, isValid);
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDxa(
				const NOX::Abstract::Vector& yVector,
				const NOX::Abstract::Vector& zVector,
				double w,
				const NOX::Abstract::MultiVector& aVector,
				NOX::Abstract::MultiVector& result_real,
				NOX::Abstract::MultiVector& result_imag)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDCeDxa(*this, yVector, zVector, w, aVector,
		  result_real, result_imag);
}

NOX::Abstract::Group::ReturnType
LOCA::Hopf::MooreSpence::FiniteDifferenceGroup::computeDCeDxa(
				const NOX::Abstract::Vector& yVector,
				const NOX::Abstract::Vector& zVector,
				double w,
				const NOX::Abstract::MultiVector& aVector,
				const NOX::Abstract::Vector& Ce_real,
				const NOX::Abstract::Vector& Ce_imag,
				NOX::Abstract::MultiVector& result_real,
				NOX::Abstract::MultiVector& result_imag)
{
  return LOCA::MultiContinuation::FiniteDifferenceGroup::derivPtr->
    computeDCeDxa(*this, yVector, zVector, w, aVector, Ce_real, Ce_imag,
		  result_real, result_imag);
}
