// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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

ostream& LOCA::StatusTest::Wrapper::print(ostream& stream, int indent) const
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
