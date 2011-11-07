#ifndef MUELU_TRILINOSSMOOTHER_DECL_HPP
#define MUELU_TRILINOSSMOOTHER_DECL_HPP

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
    TrilinosSmoother(const Xpetra::UnderlyingLib lib, std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), LO const &overlap=0, RCP<FactoryBase> AFact = Teuchos::null);
    
    //! Destructor
    virtual ~TrilinosSmoother();
    
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //! @name Setup and Apply methods.
    //@{

    //! TrilinosSmoother cannot be turned into a smoother using Setup(). Setup() always returns a RuntimeError exception.
    void Setup(Level &currentLevel);

    //! TrilinosSmoother cannot be applied. Apply() always returns a RuntimeError exception.
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const;

    //@}

    //! When this prototype is cloned using Copy(), the clone is an Ifpack or an Ifpack2 smoother.
    RCP<SmootherPrototype> Copy() const;

    //! Convert an Ifpack2 preconditioner name to Ifpack
    // As a temporary solution.
    // See https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c5 for what I proposed to do
    static std::string Ifpack2ToIfpack1Type(std::string const & type);

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
    static Teuchos::ParameterList Ifpack2ToIfpack1Param(Teuchos::ParameterList const & ifpack2List);

    //! @name Overridden from Teuchos::Describable 
    //@{
    
    //! Return a simple one-line description of this object.
    std::string description() const;
    
    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    //@}

  private:
    //! Tpetra or Epetra?
    Xpetra::UnderlyingLib lib_;

    //! ifpack1/2-specific key phrase that denote smoother type
    std::string type_;
    
    //! parameter list that is used by Ifpack/Ifpack2 internally
    Teuchos::ParameterList paramList_;

    //! overlap when using the smoother in additive Schwarz mode
    LO overlap_;

    //! A Factory
    RCP<FactoryBase> AFact_;

  }; // class TrilinosSmoother

} // namespace MueLu

#define MUELU_TRILINOS_SMOOTHER_SHORT
#endif // MUELU_TRILINOSSMOOTHER_DECL_HPP
