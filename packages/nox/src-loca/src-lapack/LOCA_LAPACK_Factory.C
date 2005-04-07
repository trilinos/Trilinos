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

#include "NOX_Parameter_List.H"

#include "LOCA_LAPACK_Factory.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_Eigensolver_DGGEVStrategy.H"

LOCA::LAPACK::Factory::Factory() :
  globalData(),
  sublistParser()
{
}

LOCA::LAPACK::Factory::~Factory()
{
}

void
LOCA::LAPACK::Factory::init(
		   const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data)
{
  globalData = global_data;
  sublistParser = Teuchos::rcp(new LOCA::Parameter::SublistParser(globalData));
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Factory::reset(
	const Teuchos::RefCountPtr<NOX::Parameter::List>& topLevelParams)
{
  // Parse sublists
  sublistParser->parseSublists(topLevelParams);

  return NOX::Abstract::Group::Ok;
}

bool
LOCA::LAPACK::Factory::createEigensolverStrategy(
	 Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy>& strategy)
{
  // Get parameter lists
  Teuchos::RefCountPtr<NOX::Parameter::List> eigenParams = 
    sublistParser->getSublist("Eigensolver");
  Teuchos::RefCountPtr<NOX::Parameter::List> solverParams = 
    sublistParser->getSublist("Linear Solver");

  // Get name of strategy
  string name = eigenParams->getParameter("Method", "Default");

  // Instantiate DGGEV strategy if requested
  if (name == "DGGEV") {
    strategy = 
      Teuchos::rcp(new LOCA::Eigensolver::DGGEVStrategy(globalData,
							eigenParams,
							solverParams));
    return true;
  }
  else
    return false;
}

