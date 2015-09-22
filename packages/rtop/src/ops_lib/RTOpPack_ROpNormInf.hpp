// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
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
