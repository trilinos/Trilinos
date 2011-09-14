#ifndef MUELU_TRILINOS_SMOOTHER_HPP
#define MUELU_TRILINOS_SMOOTHER_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherPrototype.hpp"

#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Ifpack2Smoother.hpp"

// Note: TrilinosSmoother is a SmootherPrototype that cannot be turned into a smoother using Setup().
//       When this prototype is cloned using Copy(), the clone is an Ifpack or an Ifpack2 smoother.
//       The clone can be used as a smoother after calling Setup().

namespace MueLu {

  /*!
    @class TrilinosSmoother
    @brief Class that encapsulates external library smoothers. Autoselection of Ifpack or Ifpack2 according to the underlying linear algebra library.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class TrilinosSmoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {
    
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    //! @brief Constructor
    TrilinosSmoother(const Xpetra::UnderlyingLib lib, std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), LO const &overlap=0)
      : lib_(lib), type_(type), paramList_(paramList), overlap_(overlap)
    { 
      TEST_FOR_EXCEPTION(lib_ != Xpetra::UseTpetra && lib_ != Xpetra::UseEpetra, Exceptions::RuntimeError, "lib_ != UseTpetra && lib_ != UseEpetra");
      TEST_FOR_EXCEPTION(overlap_ < 0, Exceptions::RuntimeError, "overlap_ < 0");

      if (lib_ == Xpetra::UseEpetra) {
        type_ = TrilinosSmoother::Ifpack2ToIfpack1Type(type_);
        paramList_ = TrilinosSmoother::Ifpack2ToIfpack1Param(paramList_);
      }
    }
    
    //! Destructor
    virtual ~TrilinosSmoother() {}
    
    //@}

    //! @name Setup and Apply methods.
    //@{

    //! TrilinosSmoother cannot be turned into a smoother using Setup(). Setup() always returns a RuntimeError exception.
    void Setup(Level &level) {
      TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::TrilinosSmoother::Setup(): TrilinosSmoother objects are only prototypes and TrilinosSmoother::Setup() cannot be called. Use Copy() to create an Ifpack or Ifpack2 smoother.");
    }

    //! TrilinosSmoother cannot be applied. Apply() always returns a RuntimeError exception.
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const {
      TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::TrilinosSmoother::Apply(): Setup() has not been called");
    }

    //@}

    //! When this prototype is cloned using Copy(), the clone is an Ifpack or an Ifpack2 smoother.
    RCP<SmootherPrototype> Copy() const {
      if (lib_ == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_IFPACK2
        return rcp( new Ifpack2Smoother(type_, paramList_, overlap_) );
#else
        TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external library availables for preconditionning Epetra matrices. Compile MueLu with Ifpack.");
#endif
      } else if (lib_ == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_IFPACK
        return rcp( new IfpackSmoother(type_, paramList_, overlap_) );
#else
        TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external library availables for preconditionning Tpetra matrices. Compile MueLu with Ifpack2.");
#endif
      } else {
        TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "lib_ != UseTpetra && lib_ != UseEpetra");
      }

      TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "");
    }

    //! Convert an Ifpack2 preconditioner name to Ifpack
    // As a temporary solution.
    // See https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c5 for what I proposed to do
    static std::string Ifpack2ToIfpack1Type(std::string const & type) {
      if (type == "RELAXATION") { return "point relaxation stand-alone"; }
      if (type == "CHEBYSHEV")  { return "Chebyshev"; }
    
      TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Cannot convert Ifpack2 preconditioner name to Ifpack: unkown type: " + type);
    }

    //! Convert an Ifpack2 parameter list to Ifpack
    // As a temporary solution.
    // See https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c5 for what I proposed to do
    /* TODO:
      ifpackList.set("type", "Chebyshev");
      ifpackList.set("chebyshev: degree", (int) 1);
      ifpackList.set("chebyshev: max eigenvalue", (double) 2.0);
      ifpackList.set("chebyshev: min eigenvalue", (double) 1.0);
      ifpackList.set("chebyshev: zero starting solution", false);
    */
    static Teuchos::ParameterList Ifpack2ToIfpack1Param(Teuchos::ParameterList const & ifpack2List) {
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

  private:
    //! Tpetra or Epetra?
    Xpetra::UnderlyingLib lib_;

    //! ifpack1/2-specific key phrase that denote smoother type (same as SmootherBase::Type_)
    std::string type_;
    
    //! parameter list that is used by Ifpack/Ifpack2 internally
    Teuchos::ParameterList paramList_;

    //! overlap when using the smoother in additive Schwarz mode
    LO overlap_;

  }; // class TrilinosSmoother

} // namespace MueLu

#define MUELU_TRILINOS_SMOOTHER_SHORT
#endif // MUELU_TRILINOS_SMOOTHER_HPP
