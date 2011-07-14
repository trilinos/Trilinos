#ifndef MUELU_AMESOS2_SMOOTHER_HPP
#define MUELU_AMESOS2_SMOOTHER_HPP

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities.hpp"

#ifdef HAVE_MUELU_AMESOS2
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

namespace MueLu {

template <class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class Level;

/*!
  @class Amesos2Smoother
  @brief Class that encapsulates Amesos2 direct solvers.

  This class creates an Amesos2 preconditioner factory.  The factory is capable of generating direct solvers
  based on the type and ParameterList passed into the constructor.  See the constructor for more information.
*/

  template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class Amesos2Smoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

#include "MueLu_UseShortNames.hpp"
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>   Tpetra_CrsMatrix;
    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> Tpetra_MultiVector;

  private:

    typedef Tpetra_CrsMatrix MAT; //TODO to remove
    typedef Tpetra_MultiVector MV;

    //! amesos2-specific key phrase that denote smoother type (same as SmootherBase::Type_)
    std::string amesos2Type_;
    //! pointer to Amesos2 solver object
    RCP<Amesos::Solver<MAT,MV> > prec_;
    //! matrix operator 
    Teuchos::RCP<Operator> A_;
    //! parameter list that is used by Amesos2 internally
    Teuchos::ParameterList list_;

// OLD AMESOS2
//     //! Problem that Amesos2 uses internally.
//     RCP<MultiVector> X_;
//     RCP<MultiVector> B_;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor

        Creates a MueLu interface to the direct solvers in the Amesos2 package.

    */
    Amesos2Smoother(std::string const & type, Teuchos::ParameterList const & list)
      : amesos2Type_(type), list_(list), out_(this->getOStream())
    {
      SmootherBase::SetType(type);
      SmootherPrototype::IsSetup(false);
    }

    //! Destructor
    virtual ~Amesos2Smoother() {}
    //@}

    //! @name Set/Get methods
    //@{

    //! @brief This has no effect and will throw an error.
    void SetNIts(LO const &nIts) {
      throw(Exceptions::RuntimeError("Only one iteration of Amesos2 solve is supported."));
    }

    //! @brief Returns 1.
    LO GetNIts() {
      return 1;
    }
    //@}

    //! @name Setup and Apply methods.
    //@{

    /*! @brief Set up the direct solver.

       This creates the underlying Amesos2 solver object according to the parameter list options passed into the
       Amesos2Smoother constructor.  This includes doing a numeric factorization of the matrix.
    */
    void Setup(Level &level) {
      Teuchos::OSTab tab(out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "Amesos2Smoother::Setup()" << std::endl;
      A_ = level.GetA();

// OLD AMESOS2     
//       X_ = MultiVectorFactory::Build(level.GetR()->getRangeMap(),1);
//       B_ = MultiVectorFactory::Build(level.GetR()->getRangeMap(),1);
//
      RCP<Tpetra_CrsMatrix>   tA = Utils::Op2NonConstTpetraCrs(A_);
//       RCP<Tpetra_MultiVector> tX = Utils::MV2NonConstTpetraMV(X_);
//       RCP<Tpetra_MultiVector> tB = Utils::MV2NonConstTpetraMV(B_);
//      prec_ = Amesos::create<MAT,MV>(amesos2Type_, tA, tX, tB);
  
      prec_ = Amesos::create<MAT,MV>(amesos2Type_, tA);

      if (prec_ == Teuchos::null) {
        std::string msg = "Amesos::create returns Teuchos::null";
        throw(Exceptions::RuntimeError(msg));
      }
      //TODO      prec_->setParameters(list_);

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
    void Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero=false)
    {
      if (!SmootherPrototype::IsSetup()) //TODO: use TEST_FOR_EXCEPTION
        throw(Exceptions::RuntimeError("Setup has not been called"));

// OLD AMESOS2
//       {
//         ArrayRCP<const SC> Bdata = B.getData(0);
//         ArrayRCP<SC> AmesosBdata = B_->getDataNonConst(0);
        
// //         if (*B.getMap() != *B_->getMap())
// //           throw(Exceptions::RuntimeError("Error Map Apply Amesos2Smoother"));
        
        
//         LO n = B.getMap()->getNodeNumElements();
//         for (LO i=0; i<n; i++)
//           AmesosBdata[i] = Bdata[i];
//       }

// NEW AMESOS2
      RCP<Tpetra_MultiVector> tX = Utils::MV2NonConstTpetraMV2(X);
      MultiVector & BNonC = const_cast<MultiVector&>(B);
      RCP<Tpetra_MultiVector> tB = Utils::MV2NonConstTpetraMV2(BNonC);
      prec_->setX(tX);
      prec_->setB(tB);

      //
      prec_->solve();
      //

// OLD AMESOS2
//       {
//         ArrayRCP<SC> Xdata = X.getDataNonConst(0);
//         ArrayRCP<const SC> AmesosXdata = X_->getData(0);
        
// //         if (*X.getMap() != *X_->getMap())
// //           throw(Exceptions::RuntimeError("Error Map Apply Amesos2Smoother"));
        
//         LO n = X.getMap()->getNodeNumElements();
//         for (LO i=0; i<n; i++)
//           Xdata[i] = AmesosXdata[i];
//       }
    }
    //@}

    //! @name Utilities.
    //@{

    void Print(std::string prefix) const {
      throw(Exceptions::NotImplemented("Amesos2Smoother::Print is not implemented"));
    }

    RCP<SmootherPrototype> Copy() const
    {
      return rcp(new Amesos2Smoother(*this) );
    }

    // FIXME: see Amesos2
    void CopyParameters(RCP<SmootherPrototype> source)
    {
      RCP<Amesos2Smoother> amesosSmoo = rcp_dynamic_cast<Amesos2Smoother>(source);
      //TODO check if dynamic cast fails
      amesos2Type_ = amesosSmoo->amesos2Type_;
      list_ = amesosSmoo->list_;
    }
    //@}

  }; //class Amesos2Smoother

} //namespace MueLu

#define MUELU_AMESOS2_SMOOTHER_SHORT

#endif //ifdef HAVE_MUELU_AMESOS2

#endif //ifndef MUELU_AMESOS2_SMOOTHER_HPP
