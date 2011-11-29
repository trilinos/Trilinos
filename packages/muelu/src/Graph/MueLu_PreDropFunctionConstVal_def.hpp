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
    if((Scalar)abs(vals[k]) > threshold_ || grow == gcid ) {
      //std::cout << " keep value" << std::endl;
      return false; // keep values
    }
    //std::cout << " >>>>>>>>>>>> DROP value" << std::endl;
    return true;    // values too small -> drop them
  }

}

#endif // MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP
