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

#include "LOCA_Extended_AbstractGroup.H"
#include "LOCA_Continuation_AbstractGroup.H"

const LOCA::Continuation::AbstractGroup&
LOCA::Extended::AbstractGroup::getBaseLevelUnderlyingGroup() const
{
  // First get the underlying group
  const LOCA::Continuation::AbstractGroup& ulg = getUnderlyingGroup();

  // Cast underlying group to an extended group
  const LOCA::Extended::AbstractGroup* ulgPtr = 
    dynamic_cast<const LOCA::Extended::AbstractGroup*>(&ulg);

  if (ulgPtr == NULL) {
    // Underlying group is not extended, therefore return it
    return ulg;       
  }

  else {
    // Underlying group is extended, therefore return its baselevel group
    return ulgPtr->getBaseLevelUnderlyingGroup(); 
  }

}

LOCA::Continuation::AbstractGroup&
LOCA::Extended::AbstractGroup::getBaseLevelUnderlyingGroup()
{
  // First get the underlying group
  LOCA::Continuation::AbstractGroup& ulg = getUnderlyingGroup();

  // Cast underlying group to an extended group
  LOCA::Extended::AbstractGroup* ulgPtr = 
    dynamic_cast<LOCA::Extended::AbstractGroup*>(&ulg);

  if (ulgPtr == NULL) {
    // Underlying group is not extended, therefore return it
    return ulg;       
  }

  else {
    // Underlying group is extended, therefore return its baselevel group
    return ulgPtr->getBaseLevelUnderlyingGroup(); 
  }

}
