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

#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

#include "LOCA_Bifurcation_Factory.H"
#include "LOCA_TurningPoint_MooreSpence_ExtendedGroup.H"
#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_Pitchfork_MooreSpence_ExtendedGroup.H"
#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
#include "LOCA_Pitchfork_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"

LOCA::Bifurcation::Factory::Factory(
	        const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data) : 
  globalData(global_data)
{
}

LOCA::Bifurcation::Factory::~Factory()
{
}

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
LOCA::Bifurcation::Factory::create(
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& bifurcationParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp)
{
  string methodName = "LOCA::Bifurcation::Factory::create()";
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> strategy;

  // Get name of strategy
  const string& name = strategyName(*bifurcationParams);

  if (name == "None")
    strategy = grp;

  else if (name == "Turning Point:  Moore-Spence") {

    // Cast group to MooreSpence group
    Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::AbstractGroup> msg = 
      Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MooreSpence::AbstractGroup>(grp);
    if (msg.get() == NULL)
      globalData->locaErrorCheck->throwError(
		    methodName,
		    string("Underlying group must be derived from ") + 
		    string("LOCA::TurningPoint::MooreSpence::AbstractGroup ") +
		    string("for Moore-Spence turning point continuation!"));

    strategy = 
      Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedGroup(
							   globalData,
							   topParams,
							   bifurcationParams,
							   msg));
  }
  else if (name == "Turning Point:  Minimally Augmented") {

    // Cast group to MinimallyAugmented group
    Teuchos::RefCountPtr<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup> mag = 
      Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>(grp);
    if (mag.get() == NULL)
      globalData->locaErrorCheck->throwError(
	    methodName,
	    string("Underlying group must be derived from ") + 
	    string("LOCA::TurningPoint::MinimallyAugmented::AbstractGroup ") +
	    string("for minimally augmented turning point continuation!"));

    strategy = 
      Teuchos::rcp(new LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup(
							   globalData,
							   topParams,
							   bifurcationParams,
							   mag));
  }
  else if (name == "Pitchfork:  Moore-Spence") {

    // Cast group to MooreSpence group
    Teuchos::RefCountPtr<LOCA::Pitchfork::MooreSpence::AbstractGroup> msg = 
      Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MooreSpence::AbstractGroup>(grp);
    if (msg.get() == NULL)
      globalData->locaErrorCheck->throwError(
		    methodName,
		    string("Underlying group must be derived from ") + 
		    string("LOCA::Pitchfork::MooreSpence::AbstractGroup ") +
		    string("for Moore-Spence pitchfork continuation!"));

    strategy = 
      Teuchos::rcp(new LOCA::Pitchfork::MooreSpence::ExtendedGroup(
							   globalData,
							   topParams,
							   bifurcationParams,
							   msg));
  }
  else if (name == "Pitchfork:  Minimally Augmented") {

    // Cast group to MinimallyAugmented group
    Teuchos::RefCountPtr<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup> mag = 
      Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup>(grp);
    if (mag.get() == NULL)
      globalData->locaErrorCheck->throwError(
	    methodName,
	    string("Underlying group must be derived from ") + 
	    string("LOCA::Pitchfork::MinimallyAugmented::AbstractGroup ") +
	    string("for minimally augmented pitchfork continuation!"));

    strategy = 
      Teuchos::rcp(new LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup(
							   globalData,
							   topParams,
							   bifurcationParams,
							   mag));
  }
  else if (name == "Hopf:  Moore-Spence") {

    // Cast group to MooreSpence group
    Teuchos::RefCountPtr<LOCA::Hopf::MooreSpence::AbstractGroup> msg = 
      Teuchos::rcp_dynamic_cast<LOCA::Hopf::MooreSpence::AbstractGroup>(grp);
    if (msg.get() == NULL)
      globalData->locaErrorCheck->throwError(
		    methodName,
		    string("Underlying group must be derived from ") + 
		    string("LOCA::Hopf::MooreSpence::AbstractGroup ") +
		    string("for Moore-Spence Hopf continuation!"));

    strategy = 
      Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedGroup(
							     globalData,
							     topParams,
							     bifurcationParams,
							     msg));
  }
  else if (name == "Hopf:  Minimally Augmented") {

    // Cast group to MinimallyAugmented group
    Teuchos::RefCountPtr<LOCA::Hopf::MinimallyAugmented::AbstractGroup> mag = 
      Teuchos::rcp_dynamic_cast<LOCA::Hopf::MinimallyAugmented::AbstractGroup>(grp);
    if (mag.get() == NULL)
      globalData->locaErrorCheck->throwError(
	    methodName,
	    string("Underlying group must be derived from ") + 
	    string("LOCA::Hopf::MinimallyAugmented::AbstractGroup ") +
	    string("for minimally augmented Hopf continuation!"));

    strategy = 
      Teuchos::rcp(new LOCA::Hopf::MinimallyAugmented::ExtendedGroup(
							   globalData,
							   topParams,
							   bifurcationParams,
							   mag));
  }
  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    string userDefinedName = bifurcationParams->get(
							 "User-Defined Name",
							 "???");
    if ((*bifurcationParams).INVALID_TEMPLATE_QUALIFIER
	isType< Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> >(userDefinedName))
      strategy = (*bifurcationParams).INVALID_TEMPLATE_QUALIFIER
	get< Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
				       methodName,
				       "Cannot find user-defined strategy: " + 
				       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
				      methodName,
				      "Invalid bifurcation method: " + 
				      name);

  return strategy;
}

string
LOCA::Bifurcation::Factory::strategyName(
				Teuchos::ParameterList& bifurcationParams) const
{
  // Get bifurcation type
  string bif_type =  bifurcationParams.get("Type", "None");

  // Get the formulation
  if (bif_type == "Turning Point" || bif_type == "Pitchfork" || 
      bif_type == "Hopf") {
    string formulation = 
      bifurcationParams.get("Formulation", "Moore-Spence");
    bif_type += ":  " + formulation;
  }

  return bif_type;
}
