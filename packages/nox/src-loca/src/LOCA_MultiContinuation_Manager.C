// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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

#include "NOX_Parameter_List.H"
#include "LOCA_MultiContinuation_Manager.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_ArcLengthGroup.H"
#include "LOCA_Utils.H"

LOCA::MultiContinuation::Manager::Manager(NOX::Parameter::List& p) :
  method(),
  conParamID(),
  paramsPtr(NULL)
{
  reset(p);
}

LOCA::MultiContinuation::Manager::~Manager()
{
}

NOX::Abstract::Group::ReturnType 
LOCA::MultiContinuation::Manager::reset(NOX::Parameter::List& p) 
{
  method = p.getParameter("Continuation Method", "Natural");
  conParamID = p.getParameter("Continuation Parameter", "???");
  paramsPtr = &p;
  return NOX::Abstract::Group::Ok;
}

LOCA::MultiContinuation::ExtendedGroup* 
LOCA::MultiContinuation::Manager::createContinuationGroup(
			         LOCA::MultiContinuation::AbstractGroup& grp) 
{
  if (method == "Natural") 
    return new LOCA::MultiContinuation::NaturalGroup(grp, conParamID, *paramsPtr);
  else if (method == "Arc Length")
    return new LOCA::MultiContinuation::ArcLengthGroup(grp, conParamID, *paramsPtr);
  else {
    if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
      cout << "LOCA::Continuation::Manager::createContinuationGroup() "
	   << "- invalid choice (" << method 
	   << ") for continuation method " << endl;
    }
    throw "LOCA Error";
    return NULL;
  }
}

const string&
LOCA::MultiContinuation::Manager::getMethod() const 
{
  return method;
}

const string&
LOCA::MultiContinuation::Manager::getConParamID() const 
{
  return conParamID;
}
