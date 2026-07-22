// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_SUM_HPP
#define RTOPPACK_ROP_SUM_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {


/** \brief Class ROpSum: <tt>result = sum( v0[i], i=0...n-1 )</tt> */
RTOP_ROP_1_REDUCT_SCALAR( ROpSum,
  Scalar, // Reduction object type
  REDUCT_TYPE_SUM // Reduction object reduction
  )
{
  reduct += v0;
}

} // namespace RTOpPack


#endif // RTOPPACK_ROP_SUM_HPP
