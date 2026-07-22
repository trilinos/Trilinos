// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_COUNT_NAN_INF_HPP
#define RTOPPACK_ROP_COUNT_NAN_INF_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Reduction operator that counts the number of entries that are NaN
 * or Inf.
 */
RTOP_ROP_1_REDUCT_SCALAR( ROpCountNanInf,
  index_type, // Reduction object type
  REDUCT_TYPE_SUM // Reduction object reduction
  )
{
  reduct += ( ScalarTraits<Scalar>::isnaninf(v0) ? 1 : 0 );
}


} // namespace RTOpPack


#endif // RTOPPACK_ROP_COUNT_NAN_INF_HPP
