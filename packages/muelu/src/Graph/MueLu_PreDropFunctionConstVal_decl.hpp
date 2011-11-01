#ifndef MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP
#define MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

/*
 * MueLu_PreDrop.hpp
 *
 *  Created on: Oct 26, 2011
 *      Author: agerste
 */

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

    Scalar threshold_;

  public:
    // Constructor
    explicit PreDropFunctionConstVal(const Scalar threshold);

    // Destructor
    ~PreDropFunctionConstVal() ;

    // Drop
    RCP<Graph> Drop(RCP<Operator> A) ;
  };

}

#define MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP
