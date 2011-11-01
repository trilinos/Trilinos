//TMP
#include "MueLu_PreDropFunctionBaseClass_def.hpp"
#define MUELU_PREDROPFUNCTIONBASECLASS_DECL_HPP

#ifndef MUELU_PREDROPFUNCTIONBASECLASS_DECL_HPP
#define MUELU_PREDROPFUNCTIONBASECLASS_DECL_HPP

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
    PreDropFunctionBaseClass() ;
    /// Destructor
    virtual ~PreDropFunctionBaseClass() ;

    /// Drop
    virtual RCP<Graph> Drop(RCP<Operator> A)=0;
  };
}

#define MUELU_PREDROPFUNCTIONBASECLASS_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_PREDROPFUNCTIONBASECLASS_DECL_HPP
