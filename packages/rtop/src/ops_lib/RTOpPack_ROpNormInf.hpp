// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_ROP_NORMINF_HPP
#define RTOPPACK_ROP_NORMINF_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


//RTOP_ROP_1_REDUCT_SCALAR_MAG( ROpNormInf )
//{
//  typedef ScalarTraits<Scalar> ST;
//  reduct = std::max(reduct, ST::magnitude(v0));
//}
//
//
//RTOP_ROP_1_REDUCT_SCALAR_MAG_REDUCT_OBJ_REDUCTION( ROpNormInf )
//{
//  inout_reduct = std::max(inout_reduct, in_reduct);
//}


/** \brief Norm Inf: <tt>result = max(|x[i]|, i=0...n-1). */
RTOP_ROP_1_REDUCT_SCALAR(
  ROpNormInf, // Name of the RTOp subclass
  typename ScalarTraits<Scalar>::magnitudeType, // Reduction object type
  REDUCT_TYPE_MAX // Basic reduction of reduction objects
  )
{
  typedef ScalarTraits<Scalar> ST;
  reduct = std::max(reduct, ST::magnitude(v0));
}


} // namespace RTOpPack


#endif // RTOPPACK_ROP_NORMINF_HPP
