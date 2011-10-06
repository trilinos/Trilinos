#ifndef MUELU_DIRECT_SOLVER_HPP
#define MUELU_DIRECT_SOLVER_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherPrototype.hpp"

#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Amesos2Smoother.hpp"

// Note: DirectSolver is a SmootherPrototype that cannot be turned into a smoother using Setup().
//       When this prototype is cloned using Copy(), the clone is an Amesos or an Amesos2 smoother.
//       The clone can be used as a smoother after calling Setup().

namespace MueLu {

  /*!
    @class DirectSolver
    @brief Class that encapsulates direct solvers. Autoselection of AmesosSmoother or Amesos2Smoother according to the compile time configuration of Trilinos
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class DirectSolver : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {
    
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    //! @brief Constructor
    //! Note: only parameters shared by Amesos and Amesos2 should be used for type and paramList (example: type= "Klu", "Superlu", paramList = <empty>) .
    DirectSolver(const Xpetra::UnderlyingLib lib, std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), RCP<FactoryBase> AFact = Teuchos::null)
      : lib_(lib), type_(type), paramList_(paramList), AFact_(AFact)
    { 
      TEST_FOR_EXCEPTION(lib_ != Xpetra::UseTpetra && lib_ != Xpetra::UseEpetra, Exceptions::RuntimeError, "lib_ != UseTpetra && lib_ != UseEpetra");
    }
    
    //! Destructor
    virtual ~DirectSolver() {}
    
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const {
        currentLevel.DeclareInput("A", AFact_.get());
    }

    //@}

    //! @name Setup and Apply methods.
    //@{

    //! DirectSolver cannot be turned into a smoother using Setup(). Setup() always returns a RuntimeError exception.
    void Setup(Level &currentLevel) {
      TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::DirectSolver::Setup(): DirectSolver objects are only prototypes and DirectSolver::Setup() cannot be called. Use Copy() to create an Amesos or Amesos2 smoother.");
    }

    //! DirectSolver cannot be applied. Apply() always returns a RuntimeError exception.
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const
    {
      TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::DirectSolver::Apply(): Setup() has not been called");
    }

    //@}

    //! When this prototype is cloned using Copy(), the clone is an Amesos or an Amesos2 smoother.
    RCP<SmootherPrototype> Copy() const {
      if (lib_ == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_AMESOS2
        return rcp( new Amesos2Smoother(type_, paramList_, AFact_) );
#else
        TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external direct solver library availables for Tpetra matrices. Compile MueLu with Amesos2");
#endif
      } else if (lib_ == Xpetra::UseEpetra) {
        //#if defined(HAVE_MUELU_AMESOS2)
        // return rcp( new Amesos2Smoother(type_, paramList_, AFact_) ); TODO: Amesos2 can also handle Epetra matrices but Amesos2Smoother can't for the moment.
        //#elif 

#if defined(HAVE_MUELU_AMESOS)
        return MueLu::GetAmesosSmoother<SC,LO,GO,NO,LMO>(type_, paramList_, AFact_);
#else
        TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "No external direct solver library availables for Epetra matrices. Compile MueLu with Amesos"); // add Amesos2 to the msg when done.
        return Teuchos::null;
#endif

      } else {
        TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "lib_ != UseTpetra && lib_ != UseEpetra");
      }

      TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "");
    }

    //! @name Overridden from Teuchos::Describable 
    //@{
    
    //! Return a simple one-line description of this object.
    std::string description() const {
      std::ostringstream out;
      out << SmootherPrototype::description();
      out << "{lib = " << toString(lib_) << ", type = " << type_ << "}";
      return out.str();
    }
    
    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
      MUELU_DESCRIBE;

      if (verbLevel & Parameters0) {
        out0 << "Prec. type: " << type_ << endl;
      }
      
      if (verbLevel & Parameters1) { 
        out0 << "Linear Algebra: " << toString(lib_) << endl;
        out0 << "Parameter list: " << endl; { Teuchos::OSTab tab2(out); out << paramList_; }
      }
      
      if (verbLevel & Debug) {
        out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << endl;
      }
    }

    //@}

  private:
    //! Tpetra or Epetra?
    Xpetra::UnderlyingLib lib_;

    //! amesos1/2-specific key phrase that denote smoother type
    std::string type_;
    
    //! parameter list that is used by Amesos internally
    Teuchos::ParameterList paramList_;

    //! A Factory
    RCP<FactoryBase> AFact_;

  }; // class DirectSolver

} // namespace MueLu

#define MUELU_DIRECT_SOLVER_SHORT
#endif // MUELU_DIRECT_SOLVER_HPP
