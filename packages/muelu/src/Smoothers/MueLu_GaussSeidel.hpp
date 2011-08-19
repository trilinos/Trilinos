#ifndef MUELU_GAUSSSEIDEL_HPP
#define MUELU_GAUSSSEIDEL_HPP

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {

template
<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class GaussSeidel : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

  private:

    //! sweeps
    LO nIts_;
    //! relaxation parameter
    SC omega_;
    RCP<Operator> A_;

  public:

  //! @name Constructors/Destructors
  //@{
  GaussSeidel(LO sweeps=1, SC omega=1.0) : nIts_(sweeps), omega_(omega) {
    // TODO SmootherBase::SetType(type);
    SmootherPrototype::IsSetup(false);
  }

  virtual ~GaussSeidel() {}
  //@}

  //! @name Setup and apply methods.
  //@{

  //! Set up the smoother. (Right now, just grab A from the Level.)
  void Setup(RCP<Level> const level)
  {
    A_ = level.Get< RCP<Operator> >("A");
    SmootherPrototype::IsSetup(true);
  }

  /*! Solve A*x=b approximately with the smoother.

      FIXME right now, InitialGuessIsZero is ignored.

      @param x  unknown vector
      @param b  right-hand side
      @param InitialGuessIsZero if true, indicates that x is zero, and that some flops might be avoided
  */
  void Apply(RCP<MultiVector> x, RCP<MultiVector> const b, bool InitialGuessIsZero=false)
  {
     if (InitialGuessIsZero)
       throw(Exceptions::NotImplemented("No logic for handling zero initial guesses"));

     // get matrix diagonal
     RCP<Vector> diag = VectorFactory::Build(A_->getRangeMap());
     A_->getLocalDiagCopy(*diag);

     Teuchos::ArrayRCP<const SC> bdata = b->getData(0);
     Teuchos::ArrayRCP<SC>       xdata = x->getDataNonConst(0);
     Teuchos::ArrayRCP<const SC> diagdata = diag->getData(0);

     // loop through rows
     SC sum;
     Teuchos::ArrayRCP<const LO> indices;
     Teuchos::ArrayRCP<const SC>  values;
     for (int i=0; i<A_->getNodeNumRows(); i++) {
       A_->getLocalRowView(i,indices,values);
       sum = bdata[i];
       for (int j=0; j< indices.size(); j++)
         sum -= values[j]*xdata[indices[j]];
       xdata[i] += (omega_/diagdata[i]) * sum;
     }
  } //Apply()

  //@}

  //! @name Set/Get methods.
  //@{

  void SetNIts(LO Nits) {
    nIts_ = Nits;
  }

  LO GetNIts() {
    return nIts_;
  }

  //@}

  //! @name Utilities.
  //@{

  void Print(std::string prefix) const {
    throw(Exceptions::NotImplemented("Printing not implemented yet"));
  }

  RCP<SmootherPrototype> Copy() const
  {
    return rcp(new GaussSeidel(*this) );
  }

  void CopyParameters(RCP<SmootherPrototype> source)
  {
    nIts_ = source.nIts_;
    omega_ = source.omega_;
    A_ = source.A_;
  }

  //@}

}; //class GaussSeidel

} //namespace MueLu

#define MUELU_GAUSSSEIDEL_SHORT
#endif //ifndef MUELU_GAUSSSEIDEL_HPP
