#ifndef MUELU_GAUSSSEIDELSMOOTHER_DEF_HPP
#define MUELU_GAUSSSEIDELSMOOTHER_DEF_HPP

#include <Xpetra_VectorFactory.hpp>

#include "MueLu_GaussSeidelSmoother_decl.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GaussSeidelSmoother(LO sweeps, SC omega) : nSweeps_(sweeps), omega_(omega) {
    TEUCHOS_TEST_FOR_EXCEPTION(sweeps != 1, Exceptions::NotImplemented, "MueLu::GaussSeidelSmoother(): Sweeps != 1 not implemented. Use MueLu::TrilinosSmoother instead.");
    SmootherPrototype::IsSetup(false);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  GaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~GaussSeidelSmoother() { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void GaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", NULL); //FIXME AFact_.get());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void GaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level & level) {
    A_ = level.Get< RCP<Operator> >("A", NULL); //FIXME AFact_.get());
    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void GaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &x, MultiVector const &rhs, bool const &InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(A_->getRowMap()->getComm()->getSize() != 1, Exceptions::NotImplemented, "MueLu::GaussSeidelSmoother(): This smoother is implemented only for sequential run. Use MueLu::TrilinosSmoother instead.");

    if (InitialGuessIsZero) // TODO: There is no optimization for InitialGuessIsZero = true
      x.putScalar(0.0);
      
    // get matrix diagonal
    RCP<Vector> diag = VectorFactory::Build(A_->getRangeMap());
    A_->getLocalDiagCopy(*diag);

    Teuchos::ArrayRCP<const SC> bData    = rhs.getData(0);
    Teuchos::ArrayRCP<SC>       xData    = x.getDataNonConst(0);
    Teuchos::ArrayRCP<const SC> diagData = diag->getData(0);

    // loop through rows
    SC sum;
    Teuchos::ArrayView<const LO> indices;
    Teuchos::ArrayView<const SC>  values;

    for (size_t i = 0; i < A_->getNodeNumRows(); i++) {
      A_->getLocalRowView(i, indices, values);

      sum = bData[i];
      for (int j = 0; j < indices.size(); j++)
        sum -= values[j] * xData[indices[j]];

      xData[i] += (omega_ / diagData[i]) * sum;
    }

  } // Apply ()

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > GaussSeidelSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp(new GaussSeidelSmoother(*this));
  }

} // namespace MueLu

#endif // MUELU_GAUSSSEIDELSMOOTHER_DEF_HPP
