#ifndef MUELU_GAUSSSEIDELSMOOTHER_HPP
#define MUELU_GAUSSSEIDELSMOOTHER_HPP

#include <Xpetra_Operator.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels < void, LocalOrdinal, Node>::SparseOps>
  class GaussSeidelSmoother : public SmootherPrototype <Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors/Destructors
    //@{

    GaussSeidelSmoother(LO sweeps = 1, SC omega = 1.0) : nSweeps_(sweeps), omega_(omega) {
      TEUCHOS_TEST_FOR_EXCEPTION(sweeps != 1, Exceptions::NotImplemented, "MueLu::GaussSeidelSmoother(): Sweeps != 1 not implemented. Use MueLu::TrilinosSmoother instead.");
      SmootherPrototype::IsSetup(false);
    }

    virtual ~GaussSeidelSmoother() { }

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const {
      currentLevel.DeclareInput("A", NULL); //FIXME AFact_.get());
    }

    //@}

    //! @name Setup and apply methods.
    //@{

    //! Set up the smoother. (Right now, just grab A from the Level.)
    void Setup(Level & level) {
      A_ = level.Get< RCP<Operator> >("A", NULL); //FIXME AFact_.get());
      SmootherPrototype::IsSetup(true);
    }

    /*! Solve A*x = b approximately with the smoother.
      @param x  unknown vector
      @param b  right-hand side
      @param InitialGuessIsZero if true, indicates that x is zero, and that some flops might be avoided
    */
    void Apply(MultiVector &x, MultiVector const &rhs, bool const &InitialGuessIsZero = false) const {
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
    
    //@}

    //! @name Utilities.
    //@{

    RCP <SmootherPrototype> Copy() const {
      return rcp(new GaussSeidelSmoother(*this));
    }

    //@}

  private:

    LO nSweeps_;       // < ! sweeps
    SC omega_;         // < ! relaxation parameter
    RCP <Operator> A_;

  }; //class GaussSeidelSmoother

} //namespace MueLu

#define MUELU_GAUSSSEIDELSMOOTHER_SHORT
#endif //ifndef MUELU_GAUSSSEIDELSMOOTHER_HPP
