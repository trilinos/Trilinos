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

#include "LOCA_Extended_MultiAbstractGroup.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"

Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
LOCA::Extended::MultiAbstractGroup::getBaseLevelUnderlyingGroup() const
{
  // First get the underlying group
  Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup> ulg =
    getUnderlyingGroup();

  // Cast underlying group to an extended group
  Teuchos::RCP<const LOCA::Extended::MultiAbstractGroup> ulgPtr =
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

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::Extended::MultiAbstractGroup::getBaseLevelUnderlyingGroup()
{
  // First get the underlying group
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> ulg =
    getUnderlyingGroup();

  // Cast underlying group to an extended group
  Teuchos::RCP<LOCA::Extended::MultiAbstractGroup> ulgPtr =
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
