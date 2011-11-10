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
  RCP<MueLu::Graph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Drop(RCP<Operator> A) {

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

}

#endif // MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP
