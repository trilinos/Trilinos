#ifndef MUELU_DIRECTSOLVER_DEF_HPP
#define MUELU_DIRECTSOLVER_DEF_HPP

#include <Xpetra_Utils.hpp>
#include <Xpetra_Operator.hpp>

#include "MueLu_DirectSolver_decl.hpp"

#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DirectSolver(std::string const & type, Teuchos::ParameterList const & paramList, RCP<FactoryBase> AFact)
    : type_(type), paramList_(paramList), AFact_(AFact)
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get()); // TODO: also call Amesos or Amesos2::DeclareInput?
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    //Monitor m(*this, "Setup Smoother");
    if (SmootherPrototype::IsSetup() == true) VerboseObject::GetOStream(Warnings0, 0) << "Warning: MueLu::DirectSolver::Setup(): Setup() has already been called";
    TEUCHOS_TEST_FOR_EXCEPTION(s_ != Teuchos::null, Exceptions::RuntimeError, "IsSetup() == false but s_ != Teuchos::null. This does not make sense");

    Xpetra::UnderlyingLib lib = currentLevel.Get< RCP<Operator> >("A", AFact_.get())->getRowMap()->lib();

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_AMESOS2
      s_ = rcp( new Amesos2Smoother(type_, paramList_, AFact_) );
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external direct solver library availables for Tpetra matrices. Compile MueLu with Amesos2");
#endif
    } else if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_AMESOS
      s_ = GetAmesosSmoother<SC,LO,GO,NO,LMO>(type_, paramList_, AFact_);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external direct solver library availables for Epetra matrices. Compile MueLu with Amesos"); // add Amesos2 to the msg when support for Amesos2+Epetra is implemented.
#endif
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "lib != UseTpetra && lib != UseEpetra");
    }
    
    TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError, "");

    s_->Setup(currentLevel);

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");
    TEUCHOS_TEST_FOR_EXCEPTION(s_ == Teuchos::null, Exceptions::RuntimeError, "IsSetup() == true but s_ == Teuchos::null. This does not make sense");

    s_->Apply(X, B, InitialGuessIsZero);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp( new DirectSolver(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }
    
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    //TODO
    //     if (s_ != Teuchos::null) {
    //       // Teuchos::OSTab tab2(out);
    //       s_->print(out, verbLevel);
    //     }
    
    //     if (verbLevel & Debug) {
    MUELU_DESCRIBE;
    
    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }
    
    if (verbLevel & Parameters1) { 
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab3(out); out << paramList_; }
    }

    if (verbLevel & Debug) {    
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu

#define MUELU_DIRECT_SOLVER_SHORT
#endif // MUELU_DIRECTSOLVER_DEF_HPP
