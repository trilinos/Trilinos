// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_NORM1_HPP
#define RTOPPACK_ROP_NORM1_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Class ROpNorm1. */
RTOP_ROP_1_REDUCT_SCALAR( ROpNorm1,
  typename ScalarTraits<Scalar>::magnitudeType, // Reduction object type
  REDUCT_TYPE_SUM // Reduction object reduction
  )
{
  reduct += ScalarTraits<Scalar>::magnitude(v0);
}


} // namespace RTOpPack


#endif // RTOPPACK_ROP_NORM1_HPP
