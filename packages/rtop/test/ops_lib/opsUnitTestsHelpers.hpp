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


#include "RTOpPack_RTOpTHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_as.hpp"


// Size of the vectors
extern int n;
// Cushion for machine eps
extern double errorTolSlack;
// verbose
extern bool verbose;

// 2008/07/03: rabartl: Above, we are defining these in the global namespace
// but that should be fine since these are just used for a unit test program
// and should not collide with any well-written library (which would never do
// something like this).


namespace {


using Teuchos::RCP;
using Teuchos::as;
using Teuchos::tuple;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::ScalarTraits;
using Teuchos::FancyOStream;
using RTOpPack::ScalarIndex;
using RTOpPack::ReductTarget;
using RTOpPack::DefaultReductTarget;
using RTOpPack::SubVectorView;
using RTOpPack::ConstSubVectorView;
typedef RTOpPack::index_type index_type;


template<class Scalar>
SubVectorView<Scalar>
newSubVectorView(const int n, const Scalar &val)
{
  ArrayRCP<Scalar> vals = Teuchos::arcp<Scalar>(n);
  std::fill(vals.begin(), vals.end(), val);
  return SubVectorView<Scalar>(
    0, n, vals, 1);
}


template<class Scalar>
SubVectorView<Scalar>
newStridedSubVectorView(const int n, const int stride, const Scalar &val)
{
  ArrayRCP<Scalar> vals = Teuchos::arcp<Scalar>(n*stride);
  std::fill(vals.begin(), vals.end(), Teuchos::ScalarTraits<Scalar>::nan());
  for (
    typename ArrayRCP<Scalar>::iterator itr = vals.begin();
    itr != vals.end();
    itr += stride
    )
  {
    *itr = val;
  }
  return SubVectorView<Scalar>(
    0, n, vals, stride);
}


} // namespace
