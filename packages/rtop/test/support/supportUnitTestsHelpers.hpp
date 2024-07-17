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


namespace TestingSupportHelpers {


using Teuchos::null;
using Teuchos::RCP;
using Teuchos::as;
using Teuchos::tuple;
using Teuchos::outArg;
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
void print(const ConstSubVectorView<Scalar> &sv, const std::string &sv_name,
  std::ostream &out)
{
  out << sv_name << " = " << sv << "\n";
  for (int i = 0; i < sv.subDim(); ++i) {
    out << sv[i] << ":" << i << "\n";
  }
}


template<class Scalar>
SubVectorView<Scalar>
newStridedSubVectorView(const int m, const int stride, const Scalar &val)
{
  ArrayRCP<Scalar> vals = Teuchos::arcp<Scalar>(m*stride);
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
    0, m, vals, stride);
}


template<class Scalar>
SubVectorView<Scalar>
newStridedRandomSubVectorView(const int m, const int stride)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  ArrayRCP<Scalar> vals = Teuchos::arcp<Scalar>(m*stride);
  std::fill(vals.begin(), vals.end(), ST::nan());
  for (
    typename ArrayRCP<Scalar>::iterator itr = vals.begin();
    itr != vals.end();
    itr += stride
    )
  {
    *itr = ST::random();
  }
  return SubVectorView<Scalar>(
    0, m, vals, stride);
}


template<class Scalar>
SubVectorView<Scalar>
newSubVectorView(const int m, const Scalar &val)
{
  return newStridedSubVectorView(m, 1, val);
}


template<class Scalar>
class ConstSubVectorViewAsArray {
public:
  ConstSubVectorViewAsArray(const ConstSubVectorView<Scalar> &sv)
    :sv_(sv)
    {}
  int size() const { return sv_.subDim(); }
  const Scalar operator[](int i) const { return sv_(i); }
private:
  ConstSubVectorView<Scalar> sv_;
};


template<class Scalar>
const RTOpPack::RTOpT<Scalar>& rtopt( const RTOpPack::RTOpT<Scalar> &op )
{
  return op;
}


template<class Scalar>
ConstSubVectorViewAsArray<Scalar>
constSubVectorViewAsArray(const ConstSubVectorView<Scalar> &sv)
{
  return ConstSubVectorViewAsArray<Scalar>(sv);
}


template<class Scalar>
void dumpSubVectorView(
  const ConstSubVectorView<Scalar> &sv, const std::string &sv_name,
  std::ostream &out )
{
  out << sv_name << " = {";
  for (index_type k = 0; k < sv.subDim(); ++k) {
    out << sv[k];
    if (k < sv.subDim() - 1)
      out << ", ";
  }
  out << "}\n";
}


// Size of the vectors
extern int n;

// Cushion for machine eps
extern double errorTolSlack;

// verbose
extern bool verbose;


} // namespace TestingSupportHelpers 


namespace {


using TestingSupportHelpers::n;
using TestingSupportHelpers::errorTolSlack;
using TestingSupportHelpers::verbose;
using TestingSupportHelpers::newStridedSubVectorView;
using TestingSupportHelpers::newStridedRandomSubVectorView;
using TestingSupportHelpers::newSubVectorView;
using TestingSupportHelpers::ConstSubVectorViewAsArray;
using TestingSupportHelpers::rtopt;
using TestingSupportHelpers::constSubVectorViewAsArray;
using TestingSupportHelpers::dumpSubVectorView;

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::as;
using Teuchos::tuple;
using Teuchos::outArg;
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


} // namespace
