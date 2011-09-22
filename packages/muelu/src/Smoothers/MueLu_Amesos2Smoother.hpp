// TODO: PARAMETER LIST NOT TAKE INTO ACCOUNT !!!

#ifndef MUELU_AMESOS2_SMOOTHER_HPP
#define MUELU_AMESOS2_SMOOTHER_HPP

#ifdef HAVE_MUELU_AMESOS2
#include <Amesos2_config.h>
#include <Amesos2.hpp>

#include "MueLu_ConfigDefs.hpp"
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
      : type_(type), paramList_(paramList), AFact_(AFact)
    {

#if defined(HAVE_AMESOS2_SUPERLU)
      type_ = "Superlu";
#elif defined(HAVE_AMESOS2_KLU)
      type_ = "Klu";
#endif
      TEST_FOR_EXCEPTION(type_ == "", Exceptions::RuntimeError, "MueLu::Amesos2Smoother::Amesos2Smoother(): Amesos2 compiled without KLU and SuperLU. Cannot define a solver by default for this Amesos2Smoother object");

    }

    //! Destructor
    virtual ~Amesos2Smoother() {}
    //@}

    //! @name Setup and Apply methods.
    //@{

    /*! @brief Set up the direct solver.
      This creates the underlying Amesos2 solver object according to the parameter list options passed into the
      Amesos2Smoother constructor.  This includes doing a numeric factorization of the matrix.
    */
    void Setup(Level &currentLevel) {
      Monitor m(*this, "Setup Smoother");
      if (SmootherPrototype::IsSetup() == true) this->GetOStream(Warnings0, 0) << "Warning: MueLu::Amesos2Smoother::Setup(): Setup() has already been called";

      RCP<Operator> A_ = currentLevel.Get< RCP<Operator> >("A", AFact_.get());

      RCP<Tpetra_CrsMatrix> tA = Utils::Op2NonConstTpetraCrs(A_);
  
      prec_ = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>(type_, tA);
      TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "Amesos2::create returns Teuchos::null");

      //TODO      prec_->setParameters(paramList_);
      //TODO
      // int rv = prec_->numericFactorization();
      //       if (rv != 0) {
      //         std::ostringstream buf;
      //         buf << rv;
      //         std::string msg = "Amesos2_BaseSolver::NumericFactorization return value of " + buf.str(); //TODO: BaseSolver or ... ?
      //         throw(Exceptions::RuntimeError(msg));
      //       }

      SmootherPrototype::IsSetup(true);
    }

    /*! @brief Apply the direct solver.
    Solves the linear system <tt>AX=B</tt> using the constructed solver.
    @param X initial guess
    @param B right-hand side
    @param InitialGuessIsZero This option has no effect.
    */
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false) const
    {
      TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::Amesos2Smoother::Apply(): Setup() has not been called");

      RCP<Tpetra_MultiVector> tX = Utils::MV2NonConstTpetraMV2(X);
      MultiVector & BNonC = const_cast<MultiVector&>(B);
      RCP<Tpetra_MultiVector> tB = Utils::MV2NonConstTpetraMV2(BNonC);
      prec_->setX(tX);
      prec_->setB(tB);

      prec_->solve();

      prec_->setX(Teuchos::null);
      prec_->setB(Teuchos::null);
    }
    //@}

    RCP<SmootherPrototype> Copy() const {
      return rcp( new Amesos2Smoother(*this) );
    }

    //! @name Overridden from Teuchos::Describable 
    //@{
    
    //! Return a simple one-line description of this object.
    std::string description() const {
      std::ostringstream out;
      out << SmootherPrototype::description();
      out << "{type = " << type_ << "}";
      return out.str();
    }
    
    //! Print the object with some verbosity level to an FancyOStream object.
    using MueLu::Describable::describe; // overloading, not hiding
    void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
      MUELU_DESCRIBE;

      if (verbLevel & Parameters0) {
        out0 << "Prec. type: " << type_ << endl;
      }
      
      if (verbLevel & Parameters1) { 
        out0 << "Parameter list: " << endl; { Teuchos::OSTab tab2(out); out << paramList_; }
      }
      
      if (verbLevel & External) {
        if (prec_ != Teuchos::null) { Teuchos::OSTab tab2(out); out << *prec_ << std::endl; }
      }

      if (verbLevel & Debug) {
        out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << endl
             << "-" << endl
             << "RCP<prec_>: " << prec_ << std::endl;
      }
    }

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

#endif // MUELU_AMESOS2_SMOOTHER_HPP
