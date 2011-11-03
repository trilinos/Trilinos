#ifndef MUELU_TRILINOSSMOOTHER_DEF_HPP
#define MUELU_TRILINOSSMOOTHER_DEF_HPP

#include "MueLu_TrilinosSmoother_decl.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TrilinosSmoother(const Xpetra::UnderlyingLib lib, std::string const & type, Teuchos::ParameterList const & paramList, LO const &overlap, RCP<FactoryBase> AFact)
    : lib_(lib), type_(type), paramList_(paramList), overlap_(overlap), AFact_(AFact)
  { 
    TEUCHOS_TEST_FOR_EXCEPTION(lib_ != Xpetra::UseTpetra && lib_ != Xpetra::UseEpetra, Exceptions::RuntimeError, "lib_ != UseTpetra && lib_ != UseEpetra");
    TEUCHOS_TEST_FOR_EXCEPTION(overlap_ < 0, Exceptions::RuntimeError, "overlap_ < 0");
  }
    
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~TrilinosSmoother() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::TrilinosSmoother::Setup(): TrilinosSmoother objects are only prototypes and TrilinosSmoother::Setup() cannot be called. Use Copy() to create an Ifpack or Ifpack2 smoother.");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::TrilinosSmoother::Apply(): Setup() has not been called");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    if (lib_ == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_IFPACK2
      return rcp( new Ifpack2Smoother(type_, paramList_, overlap_, AFact_) );
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external library availables for preconditionning Tpetra matrices. Compile MueLu with Ifpack2.");
#endif
    } else if (lib_ == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_IFPACK
      return GetIfpackSmoother<SC,LO,GO,NO,LMO>(TrilinosSmoother::Ifpack2ToIfpack1Type(type_), TrilinosSmoother::Ifpack2ToIfpack1Param(paramList_), overlap_, AFact_);
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external library availables for preconditionning Epetra matrices. Compile MueLu with Ifpack.");
#endif
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "lib_ != UseTpetra && lib_ != UseEpetra");
    }

    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Ifpack2ToIfpack1Type(std::string const & type) {
    if (type == "RELAXATION") { return "point relaxation stand-alone"; }
    if (type == "CHEBYSHEV")  { return "Chebyshev"; }
    
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Cannot convert Ifpack2 preconditioner name to Ifpack: unkown type: " + type);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ParameterList TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Ifpack2ToIfpack1Param(Teuchos::ParameterList const & ifpack2List) {
    Teuchos::ParameterList ifpack1List = ifpack2List;

    if (ifpack2List.isParameter("relaxation: type")) {
      std::string relaxationType = ifpack2List.get<std::string>("relaxation: type");
      if (relaxationType == "Symmetric Gauss-Seidel") {
        ifpack1List.remove("relaxation: type");
        ifpack1List.set("relaxation: type", "symmetric Gauss-Seidel");
      }
    }
      
    return ifpack1List;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{lib = " << toString(lib_) << ", type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }
      
    if (verbLevel & Parameters1) { 
      out0 << "Linear Algebra: " << toString(lib_) << std::endl;
      out0 << "PrecType: " << type_ << std::endl;
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
      out0 << "Overlap: " << overlap_ << std::endl;
    }
      
    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "Epetra PrecType: " << Ifpack2ToIfpack1Type(type_) << std::endl
           << "Epetra Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << Ifpack2ToIfpack1Param(paramList_); }
    }
  }

} // namespace MueLu

#endif // MUELU_TRILINOSSMOOTHER_DEF_HPP
