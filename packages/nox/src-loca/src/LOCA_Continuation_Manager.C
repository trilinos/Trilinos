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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_Parameter_List.H"
#include "LOCA_Continuation_Manager.H"
#include "LOCA_Abstract_Group.H"
#include "LOCA_Continuation_Group.H"
#include "LOCA_Continuation_NaturalGroup.H"
#include "LOCA_Continuation_ArcLengthGroup.H"

LOCA::Continuation::Manager::Manager(NOX::Parameter::List& params) :
  method(),
  conParamID()
{
  reset(params);
}

LOCA::Continuation::Manager::~Manager()
{
}

NOX::Abstract::Group::ReturnType 
LOCA::Continuation::Manager::reset(NOX::Parameter::List& params) 
{
  method = params.getParameter("Continuation Method", "Natural");
  conParamID = params.getParameter("Continuation Parameter", "???");
  
  return NOX::Abstract::Group::Ok;
}

LOCA::Continuation::Group* 
LOCA::Continuation::Manager::createContinuationGroup(const LOCA::Abstract::Group& grp, const NOX::Parameter::List& linSolverParams) 
{
  if (method == "Natural") 
    return new LOCA::Continuation::NaturalGroup(grp, conParamID, linSolverParams);
  else if (method == "ArcLength")
    return new LOCA::Continuation::ArcLengthGroup(grp, conParamID, linSolverParams);
  else {
    cerr << "LOCA::Continuation::Manager::createContinuationGroup() - invalid choice (" << method << ") for continuation method " << endl;
      return NULL;
  }
}

const string&
LOCA::Continuation::Manager::getMethod() const 
{
  return method;
}

const string&
LOCA::Continuation::Manager::getConParamID() const 
{
  return conParamID;
}
