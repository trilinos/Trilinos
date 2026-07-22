// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_MIN_HPP
#define RTOPPACK_ROP_MIN_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Minimum element: <tt>result = min{ v0[i], i=0...n-1 }</tt>.
 */
RTOP_ROP_1_REDUCT_SCALAR_CUSTOM_DEFAULT(
  ROpMin, // Name of the RTOp subclass
  Scalar, // Reduction object type
  REDUCT_TYPE_MIN, // Basic reduction of reduction objects
  std::numeric_limits<Scalar>::max() // Custom default reduct object value
  )
{
  reduct = std::min(reduct, v0);
}


} // namespace RTOpPack


#endif // RTOPPACK_ROP_MIN_HPP
