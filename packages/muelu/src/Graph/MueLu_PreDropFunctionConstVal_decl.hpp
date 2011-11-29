#ifndef MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP
#define MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PreDropFunctionBaseClass.hpp"
#include "MueLu_PreDropFunctionConstVal_fwd.hpp"

#include "MueLu_Graph_fwd.hpp"

namespace MueLu {

  /*!
   * Example implementation for dropping values smaller then a constant threshold
   *
   */
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class PreDropFunctionConstVal : public MueLu::PreDropFunctionBaseClass<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> {
#undef MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    // Constructor
    explicit PreDropFunctionConstVal(const Scalar threshold = 0.0);

    // Destructor
    virtual ~PreDropFunctionConstVal() { }

    /*! Drop
     * @param lrow (size_t): local row index (=lrowid)
     * @param grow (GlobalOrdinal: global row id
     * @param k    (size_t): local column iterator
     * @param lcid (LocalOrdinal): local column id (=indices[k])
     * @param gcid (GlobalOrdinal): global column id
     * @param indices (ArrrayView): array of local column ids in current row (lrow)
     * @param vals (ArrayView): array of corresponding values in current row (lrow)
     * @return bool: false, if value in (lrow, lcid) shall be kept, true if it should be dropped
     */
    bool Drop(size_t lrow, GlobalOrdinal grow, size_t k, LocalOrdinal lcid, GlobalOrdinal gcid, const Teuchos::ArrayView<const LocalOrdinal> & indices, const Teuchos::ArrayView<const Scalar> & vals);

  private:

    Scalar threshold_;

  };

}

#define MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
#endif // MUELU_PREDROPFUNCTIONCONSTVAL_DECL_HPP
