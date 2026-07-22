// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_SpmdVectorSpaceUtilities.hpp"
#include "Teuchos_CommHelpers.hpp"


namespace Thyra {


Ordinal SpmdVectorSpaceUtilities::computeMapCode(
  const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim
  )
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  //
  // Here we will make a map code out of just the local sub-dimension on each
  // processor.  If each processor has the same number of local elements, then
  // the map codes will be the same and this is all you need for RTOp
  // compatibility.
  //
  const int procRank = comm.getSize ();
  Ordinal mapCode = -1;
  Ordinal localCode = localSubDim % (procRank+1) + localSubDim;
  reduceAll<Ordinal, Ordinal> (comm, REDUCE_SUM, localCode, outArg (mapCode));
  return mapCode;
}


Ordinal SpmdVectorSpaceUtilities::computeLocalOffset(
  const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim
  )
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_SUM;
  using Teuchos::scan;

  Ordinal localOffset;
  const Ordinal _localOffset = localSubDim;
  scan<Ordinal, Ordinal> (comm, REDUCE_SUM, _localOffset, outArg (localOffset));
  localOffset -= localSubDim;
  return localOffset;
}


Ordinal SpmdVectorSpaceUtilities::computeGlobalDim(
  const Teuchos::Comm<Ordinal> &comm, const Ordinal localSubDim
  )
{
  using Teuchos::outArg;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;

  Ordinal globalDim = -1;
  reduceAll<Ordinal, Ordinal> (comm, REDUCE_SUM, localSubDim, outArg (globalDim));
  return globalDim;
}


} // namespace Thyra
