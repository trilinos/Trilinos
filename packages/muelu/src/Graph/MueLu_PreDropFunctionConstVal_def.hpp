#ifndef MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP
#define MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP

/*
 * MueLu_PreDrop.hpp
 *
 *  Created on: Oct 26, 2011
 *      Author: agerste
 */

#include "Xpetra_Operator.hpp"
#include "Xpetra_CrsGraphFactory.hpp"

#include "MueLu_ConfigDefs.hpp"
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
    /// Constructor
    explicit PreDropFunctionConstVal(const Scalar threshold)
    : threshold_(threshold) {}
    /// Destructor
    ~PreDropFunctionConstVal() {}

    /// Drop
    RCP<Graph> Drop(RCP<Operator> A) {

      RCP<CrsGraph> Xgraph = Xpetra::CrsGraphFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(A->getRowMap(), 1);

      // loop over local rows
      for(size_t row=0; row<A->getNodeNumRows(); row++)
      {
        size_t nnz = A->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        A->getLocalRowView(row, indices, vals);

        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::ThresholdAFilterFactory::Build: number of nonzeros not equal to number of indices? Error.");

        Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
        size_t nNonzeros = 0;
        for(size_t i=0; i<(size_t)indices.size(); i++) {
          if((Scalar)abs(vals[i]) > threshold_ || indices[i]==(LocalOrdinal)row) {
            indout[nNonzeros] = A->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)
            nNonzeros++;
          }
        }

        indout.resize(nNonzeros);

        Xgraph->insertGlobalIndices(A->getRowMap()->getGlobalElement(row), indout.view(0,indout.size()));
      }

      Xgraph->fillComplete(A->getDomainMap(), A->getRangeMap());

      return rcp(new Graph(Xgraph, "Graph of A"));
      // dummy implementation, returns just the original graph
      return rcp(new Graph(A->getCrsGraph(), "Graph of A"));
    }
  };

}

#define MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
#endif // MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP
