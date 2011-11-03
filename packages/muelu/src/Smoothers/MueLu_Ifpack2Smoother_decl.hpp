#ifndef MUELU_IFPACK2SMOOTHER_DECL_HPP
#define MUELU_IFPACK2SMOOTHER_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_IFPACK2

#include "Ifpack2_Factory.hpp"

#include "MueLu_VerboseObject.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {


  /*!
    @class IfpackSmoother2
    @brief Class that encapsulates Ifpack2 smoothers.

    //   This class creates an Ifpack2 preconditioner factory. The factory creates a smoother based on the
    //   type and ParameterList passed into the constructor. See the constructor for more information.
    */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class Ifpack2Smoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{
    //TODO: update doc for Ifpack2. Right now, it's a copy of the doc of IfpackSmoother
    /*! @brief Constructor

    The options passed into Ifpack2Smoother are those given in the Ifpack2 user's manual.

    @param type smoother type
    @param list options for the particular smoother (e.g., fill factor or damping parameter)

    Here is how to select some of the most common smoothers.

    - Gauss-Seidel
    - <tt>type</tt> = <tt>point relaxation stand-alone</tt>
    - parameter list options
    - <tt>relaxation: type</tt> = <tt>Gauss-Seidel</tt>
    - <tt>relaxation: damping factor</tt>
    - symmetric Gauss-Seidel
    - <tt>type</tt> = <tt>point relaxation stand-alone</tt>
    - parameter list options
    - <tt>relaxation: type</tt> = <tt>symmetric Gauss-Seidel</tt>
    - <tt>relaxation: damping factor</tt>
    - Chebyshev
    - <tt>type</tt> = <tt>Chebyshev</tt>
    - parameter list options
    - <tt>chebyshev: ratio eigenvalue</tt>
    - <tt>chebyshev: min eigenvalue</tt>
    - <tt>chebyshev: max eigenvalue</tt>
    - <tt>chebyshev: degree</tt>
    - <tt>chebyshev: zero starting solution</tt> (defaults to <tt>true</tt>)
    - ILU
    - <tt>type</tt> = <tt>ILU</tt>
    - parameter list options
    - <tt>fact: level-of-fill</tt>

    See also Ifpack2_Relaxation, Ifpack2_Chebyshev, Ifpack2_ILUT.
    */
    Ifpack2Smoother(std::string const & type, Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), LO const &overlap=0, RCP<FactoryBase> AFact = Teuchos::null) //TODO: empty paramList valid for Ifpack??
    ;

    //! Destructor
    virtual ~Ifpack2Smoother() ;

    //@}

    //! @name Set/Get methods
    //@{

   //! Set smoother parameters
    void SetParameters(Teuchos::ParameterList const & paramList) ;

    //! Get smoother parameters
    Teuchos::ParameterList const & GetParameters() ;

    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const ;

    //@}

    //! @name Computational methods.
    //@{

    /*! @brief Set up the smoother.

    This creates the underlying Ifpack2 smoother object, copies any parameter list options
    supplied to the constructor to the Ifpack2 object, and computes the preconditioner.

    TODO The eigenvalue estimate should come from A_, not the Ifpack2 parameter list.
    */
    void Setup(Level &currentLevel) ;

    /*! @brief Apply the preconditioner.

    Solves the linear system <tt>AX=B</tt> using the constructed smoother.

    @param X initial guess
    @param B right-hand side
    @param InitialGuessIsZero (optional) If false, some work can be avoided. Whether this actually saves any work depends on the underlying Ifpack2 implementation.
    */
    void Apply(MultiVector& X, MultiVector const &B, bool const &InitialGuessIsZero=false) const
    ;

    //@}

    //! @name Utilities
    //@{

    RCP<SmootherPrototype> Copy() const ;

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{
    
    //! Return a simple one-line description of this object.
    std::string description() const ;
    
    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const 
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const ;

    //@}

  private:

    //! ifpack2-specific key phrase that denote smoother type
    std::string type_;

    //! parameter list that is used by Ifpack internally
    Teuchos::ParameterList paramList_;

    //! overlap when using the smoother in additive Schwarz mode
    LO overlap_;

    //! pointer to Ifpack2 preconditioner object
    RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec_;

    //! A Factory
    RCP<FactoryBase> AFact_;

  }; // class Ifpack2Smoother

} // namespace MueLu

#define MUELU_IFPACK2_SMOOTHER_SHORT
#endif //ifdef HAVE_MUELU_IFPACK2
#endif // MUELU_IFPACK2SMOOTHER_DECL_HPP
