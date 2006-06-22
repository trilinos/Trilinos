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

Index SpmdVectorSpaceUtilities::computeMapCode(
  const Teuchos::Comm<Index> &comm, const Index localSubDim
  )
{
  //
  // Here we will make a map code out of just the local
  // sub-dimension on each processor.  If each processor
  // has the same number of local elements, then the maps
  // will be the same and this is all you need for
  // RTOp compatibility unless the operations are not
  // coordinate invariant.  I will work on this issue
  // if it becomes a problem.
  //
  const int procRank = size(comm);
  Index mapCode = -1;
  Index localCode = localSubDim % (procRank+1) + localSubDim;
  reduceAll(comm,Teuchos::REDUCE_SUM,1,&localCode,&mapCode);
  return mapCode;
}

Index SpmdVectorSpaceUtilities::computeLocalOffset(
  const Teuchos::Comm<Index> &comm, const Index localSubDim
  )
{
  Index localOffset;
  Index _localOffset = localSubDim;
  scan(comm,Teuchos::REDUCE_SUM,1,&_localOffset,&localOffset);
  localOffset -= localSubDim;
  return localOffset;
}

Index SpmdVectorSpaceUtilities::computeGlobalDim(
  const Teuchos::Comm<Index> &comm, const Index localSubDim
  )
{
  Index globalDim = -1;
  reduceAll(comm,Teuchos::REDUCE_SUM,1,&localSubDim,&globalDim);
  return globalDim;
}

void SpmdVectorSpaceUtilities::broadcast(
  const Teuchos::Comm<Index> &comm, const int rootRank, Index* value
  )
{
  broadcast(comm,rootRank,1,value);
}

} // namespace Thyra
