#ifndef MUELU_DIRECTSOLVER_DEF_HPP
#define MUELU_DIRECTSOLVER_DEF_HPP

#include "MueLu_DirectSolver_decl.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DirectSolver(const Xpetra::UnderlyingLib lib, std::string const & type, Teuchos::ParameterList const & paramList, RCP<FactoryBase> AFact)
    : lib_(lib), type_(type), paramList_(paramList), AFact_(AFact)
  { 
    TEUCHOS_TEST_FOR_EXCEPTION(lib_ != Xpetra::UseTpetra && lib_ != Xpetra::UseEpetra, Exceptions::RuntimeError, "lib_ != UseTpetra && lib_ != UseEpetra");
  }
    
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~DirectSolver() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::DirectSolver::Setup(): DirectSolver objects are only prototypes and DirectSolver::Setup() cannot be called. Use Copy() to create an Amesos or Amesos2 smoother.");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::DirectSolver::Apply(): Setup() has not been called");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    if (lib_ == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_AMESOS2
      return rcp( new Amesos2Smoother(type_, paramList_, AFact_) );
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external direct solver library availables for Tpetra matrices. Compile MueLu with Amesos2");
#endif
    } else if (lib_ == Xpetra::UseEpetra) {
      //#if defined(HAVE_MUELU_AMESOS2)
      // return rcp( new Amesos2Smoother(type_, paramList_, AFact_) ); TODO: Amesos2 can also handle Epetra matrices but Amesos2Smoother can't for the moment.
      //#elif 

#if defined(HAVE_MUELU_AMESOS)
      return MueLu::GetAmesosSmoother<SC,LO,GO,NO,LMO>(type_, paramList_, AFact_);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external direct solver library availables for Epetra matrices. Compile MueLu with Amesos"); // add Amesos2 to the msg when done.
      return Teuchos::null;
#endif

    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "lib_ != UseTpetra && lib_ != UseEpetra");
    }

    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{lib = " << toString(lib_) << ", type = " << type_ << "}";
    return out.str();
  }
    
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void DirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }
      
    if (verbLevel & Parameters1) { 
      out0 << "Linear Algebra: " << toString(lib_) << std::endl;
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
    }
      
    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu

#define MUELU_DIRECT_SOLVER_SHORT
#endif // MUELU_DIRECTSOLVER_DEF_HPP
