// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_DOT_PROD_HPP
#define RTOPPACK_ROP_DOT_PROD_HPP


#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {


RTOP_ROP_2_REDUCT_SCALAR( ROpDotProd,
  Scalar, // Reduction object type
  REDUCT_TYPE_SUM // Reduction object reduction operator type
  )
{
  reduct += ScalarTraits<Scalar>::conjugate(v0)*(v1);
}


} // namespace RTOpPack


#endif // RTOPPACK_ROP_DOT_PROD_HPP
