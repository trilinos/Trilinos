/*
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
*/


#include "RTOpPack_TOpSetSubVector.hpp"

#include "supportUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( TOpSetSubVector, dense, Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

  SubVectorView<Scalar> input_sv = newSubVectorView<Scalar>(n, ST::zero());
  for (index_type k = 0; k < n; ++k)
    input_sv(k) = ST::random();
  
  RTOpPack::TOpSetSubVector<Scalar> setSubVectorOp;
  setSubVectorOp.set_sub_vec(input_sv);
  // 2008/07/23: rabartl: Above: For some reason, I have to use the default
  // constructor for this class object and then set the subvector or gcc 3.4.6
  // issues an error and thinks that setSubVectorOp is an operator function.
  // This must be a compiler bug!

  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, ST::nan());
  setSubVectorOp.apply_op( null, tuple(sv)(), null );

  if (verbose) {
    dumpSubVectorView(input_sv, "input_sv", out);
    dumpSubVectorView(sv, "sv", out);
  }

  TEST_COMPARE_ARRAYS( constSubVectorViewAsArray(input_sv),
    constSubVectorViewAsArray(sv) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( TOpSetSubVector, dense )


} // namespace
