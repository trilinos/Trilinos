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

#include "LOCA_StatusTest_Wrapper.H"
#include "LOCA_Solver_Wrapper.H"

LOCA::StatusTest::Wrapper::Wrapper(
		    const Teuchos::RCP<NOX::StatusTest::Generic>& s) :
  statusTestPtr(s)
{
}

LOCA::StatusTest::Wrapper::~Wrapper()
{
}

NOX::StatusTest::StatusType 
LOCA::StatusTest::Wrapper::checkStatus(const NOX::Solver::Generic& problem, 
				       NOX::StatusTest::CheckType checkType)
{
  LOCA::Solver::Wrapper problemWrapper(Teuchos::rcp(&problem,false));

  return statusTestPtr->checkStatus(problemWrapper, checkType);
}

NOX::StatusTest::StatusType 
LOCA::StatusTest::Wrapper::getStatus() const
{
  return statusTestPtr->getStatus();
}

std::ostream& LOCA::StatusTest::Wrapper::print(std::ostream& stream, int indent) const
{
  return statusTestPtr->print(stream, indent);
}

Teuchos::RCP<NOX::StatusTest::Generic>
LOCA::StatusTest::Wrapper::getUnderlyingStatusTest() 
{
  return statusTestPtr;
}

Teuchos::RCP<const NOX::StatusTest::Generic>
LOCA::StatusTest::Wrapper::getUnderlyingStatusTest() const
{
  return statusTestPtr;
}
