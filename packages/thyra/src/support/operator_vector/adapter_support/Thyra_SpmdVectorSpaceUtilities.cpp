// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Thyra_SpmdVectorSpaceUtilities.hpp"
#include "Teuchos_CommHelpers.hpp"


namespace Thyra {


Ordinal SpmdVectorSpaceUtilities::computeMapCode(
  const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim
  )
{
  //
  // Here we will make a map code out of just the local sub-dimension on each
  // processor.  If each processor has the same number of local elements, then
  // the map codes will be the same and this is all you need for RTOp
  // compatibility.
  //
  const int procRank = size(comm);
  Ordinal mapCode = -1;
  Ordinal localCode = localSubDim % (procRank+1) + localSubDim;
  reduceAll(comm, Teuchos::REDUCE_SUM, localCode,
    Teuchos::outArg(mapCode));
  return mapCode;
}


Ordinal SpmdVectorSpaceUtilities::computeLocalOffset(
  const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim
  )
{
  Ordinal localOffset;
  const Ordinal _localOffset = localSubDim;
  scan(comm, Teuchos::REDUCE_SUM, _localOffset, 
    Teuchos::outArg(localOffset));
  localOffset -= localSubDim;
  return localOffset;
}


Ordinal SpmdVectorSpaceUtilities::computeGlobalDim(
  const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim
  )
{
  Ordinal globalDim = -1;
  reduceAll(comm, Teuchos::REDUCE_SUM, localSubDim,
    Teuchos::outArg(globalDim));
  return globalDim;
}


} // namespace Thyra
