// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP
#define MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP

#include <Xpetra_CrsGraphFactory.hpp>

#include "MueLu_PreDropFunctionConstVal_decl.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PreDropFunctionConstVal(const Scalar threshold)
  : threshold_(threshold) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Drop(size_t /* lrow */, GlobalOrdinal grow, size_t k, LocalOrdinal /* lcid */, GlobalOrdinal gcid, const Teuchos::ArrayView<const LocalOrdinal>& /* indices */, const Teuchos::ArrayView<const Scalar>& vals) {
  if (Teuchos::ScalarTraits<Scalar>::magnitude(vals[k]) > Teuchos::ScalarTraits<Scalar>::magnitude(threshold_) || grow == gcid) {
    return false;  // keep values
  }
  return true;  // values too small -> drop them
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetThreshold() const {
  return threshold_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
  std::ostringstream out;
  out << "PreDropFunctionConstVal: threshold = " << threshold_ << std::endl;
  return out.str();
}

/*template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
void PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node>::describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
  MUELU_DESCRIBE;
  if (verbLevel & Parameters0) {
    out0 << "PreDropFunctionConstVal: threshold = " << threshold_ << std::endl;
  }
}*/

}  // namespace MueLu

#define MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
#endif  // MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP
