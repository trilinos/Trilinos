#ifndef MUELU_PREDROPFUNCTIONBASECLASS_DECL_HPP
#define MUELU_PREDROPFUNCTIONBASECLASS_DECL_HPP

#include <Xpetra_Operator_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_PreDropFunctionBaseClass_fwd.hpp"

namespace MueLu {

  /*!
   * Base class you can derive from to allow user defined dropping
   *
   */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class PreDropFunctionBaseClass : public BaseClass {
#undef MUELU_PREDROPFUNCTIONBASECLASS_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //! Destructor
    virtual ~PreDropFunctionBaseClass() { }

    //! Drop
    virtual bool Drop(size_t lrow, GlobalOrdinal grow, size_t k, LocalOrdinal lcid, GlobalOrdinal gcid, const Teuchos::ArrayView<const LocalOrdinal> & indices, const Teuchos::ArrayView<const Scalar> & vals) = 0;

  };
}

#define MUELU_PREDROPFUNCTIONBASECLASS_SHORT
#endif // MUELU_PREDROPFUNCTIONBASECLASS_DECL_HPP
