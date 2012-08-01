#ifndef MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP
#define MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP

#include <Xpetra_CrsGraphFactory.hpp>

#include "MueLu_PreDropFunctionConstVal_decl.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PreDropFunctionConstVal(const Scalar threshold)
    : threshold_(threshold) { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Drop(size_t lrow, GlobalOrdinal grow, size_t k, LocalOrdinal lcid, GlobalOrdinal gcid, const Teuchos::ArrayView<const LocalOrdinal> & indices, const Teuchos::ArrayView<const Scalar> & vals) {
    if(std::abs(vals[k]) > std::abs(threshold_) || grow == gcid ) {
      return false; // keep values
    }
    return true;    // values too small -> drop them
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << "PreDropFunctionConstVal: threshold = " << threshold_ << std::endl;
    return out.str();
  }

  /*template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;
    if (verbLevel & Parameters0) {
      out0 << "PreDropFunctionConstVal: threshold = " << threshold_ << std::endl;
    }
  }*/

}

#define MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
#endif // MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP
