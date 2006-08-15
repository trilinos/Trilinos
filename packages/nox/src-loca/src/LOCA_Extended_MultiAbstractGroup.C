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

#include "LOCA_Extended_MultiAbstractGroup.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"

Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup>
LOCA::Extended::MultiAbstractGroup::getBaseLevelUnderlyingGroup() const
{
  // First get the underlying group
  Teuchos::RefCountPtr<const LOCA::MultiContinuation::AbstractGroup> ulg = 
    getUnderlyingGroup();

  // Cast underlying group to an extended group
  Teuchos::RefCountPtr<const LOCA::Extended::MultiAbstractGroup> ulgPtr = 
    Teuchos::rcp_dynamic_cast<const LOCA::Extended::MultiAbstractGroup>(ulg);

  if (ulgPtr.get() == NULL) {
    // Underlying group is not extended, therefore return it
    return ulg;       
  }

  else {
    // Underlying group is extended, therefore return its baselevel group
    return ulgPtr->getBaseLevelUnderlyingGroup(); 
  }

}

Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
LOCA::Extended::MultiAbstractGroup::getBaseLevelUnderlyingGroup()
{
  // First get the underlying group
  Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup> ulg = 
    getUnderlyingGroup();

  // Cast underlying group to an extended group
  Teuchos::RefCountPtr<LOCA::Extended::MultiAbstractGroup> ulgPtr = 
    Teuchos::rcp_dynamic_cast<LOCA::Extended::MultiAbstractGroup>(ulg);

  if (ulgPtr.get() == NULL) {
    // Underlying group is not extended, therefore return it
    return ulg;       
  }

  else {
    // Underlying group is extended, therefore return its baselevel group
    return ulgPtr->getBaseLevelUnderlyingGroup(); 
  }

}
