// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_MAX_HPP
#define RTOPPACK_ROP_MAX_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Maximum element: <tt>result = max{ v0[i], i=0...n-1 }</tt>.
 */
RTOP_ROP_1_REDUCT_SCALAR_CUSTOM_DEFAULT(
  ROpMax, // Name of the RTOp subclass
  Scalar, // Reduction object type
  REDUCT_TYPE_MAX, // Basic reduction of reduction objects
  Teuchos::as<Scalar>(-std::numeric_limits<Scalar>::max())
  )
{
  reduct = std::max(reduct, v0);
}


} // namespace RTOpPack


#endif // RTOPPACK_ROP_MAX_HPP
