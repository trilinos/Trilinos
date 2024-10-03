// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
