#ifndef MUELU_PREDROPFUNCTIONBASECLASS_DEF_HPP
#define MUELU_PREDROPFUNCTIONBASECLASS_DEF_HPP

/*
 * MueLu_PreDrop.hpp
 *
 *  Created on: Oct 26, 2011
 *      Author: agerste
 */

#include "Xpetra_Operator.hpp"

//TMP
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Graph.hpp"

// #include "MueLu_ConfigDefs.hpp"
// #include "MueLu_SingleLevelFactoryBase.hpp"
// #include "MueLu_Level_def.hpp"
// #include "MueLu_Graph_def.hpp"

namespace MueLu {

  /*!
   * Base class you can derive from to allow user defined dropping
   *
   */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class PreDropFunctionBaseClass {

#include "MueLu_UseShortNames.hpp"

  public:
    /// Constructor
    PreDropFunctionBaseClass() {}

    /// Destructor
    virtual ~PreDropFunctionBaseClass() {}

    /// Drop
    virtual RCP<Graph> Drop(RCP<Operator> A)=0;
  };
}

#define MUELU_PREDROPFUNCTIONBASECLASS_SHORT
#endif // MUELU_PREDROPFUNCTIONBASECLASS_DEF_HPP
