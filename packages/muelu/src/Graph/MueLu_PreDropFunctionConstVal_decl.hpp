#ifndef MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP
#define MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "Xpetra_Operator.hpp"

#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_PreDropFunctionBaseClass.hpp"

namespace MueLu {

  /*!
   * Example implementation for dropping values smaller then a constant threshold
   *
   */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class PreDropFunctionConstVal : public MueLu::PreDropFunctionBaseClass<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

  public:

    // Constructor
    explicit PreDropFunctionConstVal(const Scalar threshold = 0.0);

    // Destructor
    virtual ~PreDropFunctionConstVal() { }

    // Drop
    RCP<Graph> Drop(RCP<Operator> A);

  private:

    Scalar threshold_;

  };

}

#define MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
#endif // MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP
