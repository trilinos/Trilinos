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
#include "LOCA_Bifurcation_Manager.H"
#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"
#include "LOCA_Bifurcation_TPBord_ExtendedGroup.H"
#include "LOCA_Bifurcation_TPBord_ModifiedBorderingGroup.H"
#include "LOCA_Bifurcation_TPBord_NicDayModifiedBorderingGroup.H"
#include "LOCA_Bifurcation_PitchforkBord_ExtendedGroup.H"
#include "LOCA_Bifurcation_HopfBord_AbstractGroup.H"
#include "LOCA_Bifurcation_HopfBord_ExtendedGroup.H"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"

LOCA::Bifurcation::Manager::Manager(NOX::Parameter::List& p) :
  method(),
  paramsPtr(NULL)
{
  reset(p);
}

LOCA::Bifurcation::Manager::~Manager()
{
}

NOX::Abstract::Group::ReturnType 
LOCA::Bifurcation::Manager::reset(NOX::Parameter::List& p) 
{
  if (!p.isParameter("Method")) {
    LOCA::ErrorCheck::printWarning(
	       "LOCA::Bifurcation::Manager::reset()",
	       "\"Method\"  is not set.  Defaulting to \"None\"");
  }
  method = p.getParameter("Method", "None");

  paramsPtr = &p;
  return NOX::Abstract::Group::Ok;
}

LOCA::Continuation::AbstractGroup* 
LOCA::Bifurcation::Manager::createBifurcationGroup(
			            LOCA::Continuation::AbstractGroup& grp) 
{
  if (method == "None") 
    return dynamic_cast<LOCA::Continuation::AbstractGroup*>(grp.clone(NOX::DeepCopy));
  else if (method == "Turning Point") {
    LOCA::Bifurcation::TPBord::AbstractGroup& tpGrp =
      dynamic_cast<LOCA::Bifurcation::TPBord::AbstractGroup&>(grp);
    return new LOCA::Bifurcation::TPBord::ExtendedGroup(tpGrp, *paramsPtr);
  }
  else if (method == "Modified Turning Point") {
     LOCA::Bifurcation::TPBord::AbstractGroup& tpGrp =
      dynamic_cast<LOCA::Bifurcation::TPBord::AbstractGroup&>(grp);
    return new LOCA::Bifurcation::TPBord::ModifiedBorderingGroup(tpGrp, 
								 *paramsPtr);
  }
  else if (method == "Nic-Day Modified Turning Point") {
     LOCA::Bifurcation::TPBord::AbstractGroup& tpGrp =
      dynamic_cast<LOCA::Bifurcation::TPBord::AbstractGroup&>(grp);
    return new LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup(tpGrp, 
								 *paramsPtr);
  }
  else if (method == "Pitchfork") {
     LOCA::Bifurcation::TPBord::AbstractGroup& tpGrp =
      dynamic_cast<LOCA::Bifurcation::TPBord::AbstractGroup&>(grp);
    return new LOCA::Bifurcation::PitchforkBord::ExtendedGroup(tpGrp, 
							       *paramsPtr);
  }
  else if (method == "Hopf") {
     LOCA::Bifurcation::HopfBord::AbstractGroup& hopfGrp =
      dynamic_cast<LOCA::Bifurcation::HopfBord::AbstractGroup&>(grp);
    return new LOCA::Bifurcation::HopfBord::ExtendedGroup(hopfGrp, *paramsPtr);
  }
  else {
    if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
      cout << "LOCA::Bifurcation::Manager::createBifurcationGroup() "
	   << "- invalid choice (" << method 
	   << ") for bifurcation method " << endl;
    }
    throw "LOCA Error";
    return NULL;
  }
}

const string&
LOCA::Bifurcation::Manager::getMethod() const 
{
  return method;
}
