#ifndef MUELU_AMESOS2SMOOTHER_DECL_HPP
#define MUELU_AMESOS2SMOOTHER_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

// TODO: PARAMETER LIST NOT TAKE INTO ACCOUNT !!!

#ifdef HAVE_MUELU_AMESOS2
#include <Amesos2_config.h>
#include <Amesos2.hpp>

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  /*!
    @class Amesos2Smoother
    @brief Class that encapsulates Amesos2 direct solvers.

    This class creates an Amesos2 preconditioner factory.  The factory is capable of generating direct solvers
    based on the type and ParameterList passed into the constructor.  See the constructor for more information.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class Amesos2Smoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {
    
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor
      Creates a MueLu interface to the direct solvers in the Amesos2 package.
      If you are using type=="", then either SuperLU or KLU2 are used by default.
    */
    Amesos2Smoother(std::string const & type = "", Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), RCP<FactoryBase> AFact = Teuchos::null)
    ;

    //! Destructor
    virtual ~Amesos2Smoother() ;
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const ;

    //@}

    //! @name Setup and Apply methods.
    //@{

    /*! @brief Set up the direct solver.
      This creates the underlying Amesos2 solver object according to the parameter list options passed into the
      Amesos2Smoother constructor.  This includes doing a numeric factorization of the matrix.
    */
    void Setup(Level &currentLevel) ;

    /*! @brief Apply the direct solver.
    Solves the linear system <tt>AX=B</tt> using the constructed solver.
    @param X initial guess
    @param B right-hand side
    @param InitialGuessIsZero This option has no effect.
    */
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const
    ;
    //@}

    RCP<SmootherPrototype> Copy() const ;

    //! @name Overridden from Teuchos::Describable 
    //@{
    
    //! Return a simple one-line description of this object.
    std::string description() const ;
    
    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const ;

    //@}

  private:
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Tpetra_CrsMatrix;
    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> Tpetra_MultiVector;

    //! amesos2-specific key phrase that denote smoother type
    std::string type_;

    //! parameter list that is used by Amesos2 internally
    Teuchos::ParameterList paramList_;

    //! pointer to Amesos2 solver object
    RCP<Amesos2::Solver<Tpetra_CrsMatrix,Tpetra_MultiVector> > prec_;

    //! A Factory
    RCP<FactoryBase> AFact_;

  }; // class Amesos2Smoother

} // namespace MueLu

#define MUELU_AMESOS2_SMOOTHER_SHORT

#endif // HAVE_MUELU_AMESOS2

#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_AMESOS2SMOOTHER_DECL_HPP
