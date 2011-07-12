#ifndef MUELU_IFPACK2_SMOOTHER_HPP
#define MUELU_IFPACK2_SMOOTHER_HPP

#ifdef HAVE_MUELU_IFPACK2

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities.hpp"

#include "Ifpack2_Factory.hpp"

namespace MueLu {

template <class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class Level;

/*!
  @class IfpackSmoother2
  @brief Class that encapsulates Ifpack2 smoothers.

//   This class creates an Ifpack2 preconditioner factory. The factory creates a smoother based on the
//   type and ParameterList passed into the constructor. See the constructor for more information.
*/

  template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class Ifpack2Smoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

#include "MueLu_UseShortNames.hpp"

  private:

    //! Ifpack2-specific key phrase that denote smoother type (not to be confused with SmootherBase::Type_)
    std::string ifpack2Type_;
    //! overlap when using the smoother in additive Schwarz mode
    LO overlap_;
    RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec_;
    //! matrix operator 
    Teuchos::RCP<Operator> A_;
    //! parameter list that is used by Ifpack2 internally
    Teuchos::ParameterList list_;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

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
    Ifpack2Smoother(std::string const & type, Teuchos::ParameterList & list)
      : ifpack2Type_(type), list_(list), out_(this->getOStream())
    {
      overlap_ = list.get("Overlap",(LO) 0);
      std::string label;
//       if (type == "point relaxation stand-alone")
//         label = "Ifpack2: " + list.get("relaxation: type","unknown relaxation");
//       else
      label = "Ifpack2: " + type;
      SmootherBase::SetType(label);
    }

    //! Destructor
    virtual ~Ifpack2Smoother() {}
    //@}

    //! @name Set/Get methods

    //@{

    /*! @brief Set the number of smoothing sweeps/degree.

       If the smoother is relaxation, this sets the number of sweeps.
       If the smoother is Chebyshev, this sets the polynomial degree.

       Note:  This can be called after the preconditioner is set up, i.e., after
       calling Ifpack2Smoother::Setup().
    */
    void SetNIts(LO const &nIts) {
//       if (!SmootherPrototype::IsSetup()) //FIXME precond doesn't have to be setup
//         throw(Exceptions::RuntimeError("Call Setup before setting sweeps"));
//       if (ifpackType_ == "point relaxation stand-alone") list_.set("relaxation: sweeps", nIts);
//       else if (ifpackType_ == "Chebyshev")               list_.set("chebyshev: degree", nIts);
//       prec_->SetParameters(list_);
      throw(Exceptions::NotImplemented("Not Implemented"));
    }

    /*! @brief Get the number of smoothing sweeps.

       If the smoother is relaxation, this returns the number of sweeps.
       If the smoother is Chebyshev, this returns the polynomial degree.
    */
    LO GetNIts() {
//       if (ifpackType_ == "point relaxation stand-alone")
//       {
//         if (list_.isParameter("relaxation: sweeps") == false)
//           throw(Exceptions::RuntimeError("number of iterations is not set"));
//         return list_.get("relaxation: sweeps",1);
//       } else if (ifpackType_ == "Chebyshev") {
//         if (list_.isParameter("chebyshev: degree") == false)
//           throw(Exceptions::RuntimeError("Chebyshev degree is not set"));
//         return list_.get("chebyshev: degree",1);
//       } else 
//         throw(Exceptions::RuntimeError("GetNIts: unknown smoother type"));

      throw(Exceptions::NotImplemented("Not Implemented"));
      return -1;
    }
    //@}

    //! @name Computational methods.
    //@{

    /*! @brief Set up the smoother.

        This creates the underlying Ifpack2 smoother object, copies any parameter list options
        supplied to the constructor to the Ifpack2 object, and computes the preconditioner.

        TODO The eigenvalue estimate should come from A_, not the Ifpack2 parameter list.
    */
    void Setup(Level &level) {
      Teuchos::OSTab tab(out_);
      A_ = level.GetA();

      // output information
      std::ostringstream buf; buf << level.GetLevelID();
      std::string prefix = "Smoother (level " + buf.str() + ") : ";
      LO rootRank = out_->getOutputToRootOnly();
      out_->setOutputToRootOnly(0);
      *out_ << prefix << "# global rows = " << A_->getGlobalNumRows()
            << ", estim. global nnz = " << A_->getGlobalNumEntries() << std::endl;

      RCP<const Tpetra::CrsMatrix<SC, LO, GO, NO, LMO> > tpA = Utils::Op2NonConstTpetraCrs(A_);
      prec_ = Ifpack2::Factory::create(ifpack2Type_, tpA, overlap_);
      if (ifpack2Type_ == "CHEBYSHEV") {
        Scalar maxEigenValue = list_.get("chebyshev: max eigenvalue",(Scalar)-1.0);
        if (maxEigenValue == -1.0) {
          maxEigenValue = Utils::PowerMethod(*A_,true,10,1e-4);
          list_.set("chebyshev: max eigenvalue",maxEigenValue);
        }
        *out_ << prefix << "Ifpack2 Chebyshev, degree " << list_.get("chebyshev degree",1) << std::endl;
        *out_ << prefix << "lambda_min=" << list_.get("chebyshev: min eigenvalue",-1.0)
              << ", lambda_max=" << list_.get("chebyshev: max eigenvalue",-1.0) << std::endl;
      }
      out_->setOutputToRootOnly(rootRank);
      prec_->setParameters(list_);
      prec_->initialize();
      prec_->compute();

      SmootherPrototype::IsSetup(true);
    }

    /*! @brief Apply the preconditioner.

        Solves the linear system <tt>AX=B</tt> using the constructed smoother.

        @param X initial guess
        @param B right-hand side
        @param InitialGuessIsZero (optional) If false, some work can be avoided.  Whether this actually saves any work depends on the underlying Ifpack2 implementation.
    */
    void Apply(MultiVector& X, MultiVector const &B, bool const &InitialGuessIsZero=false)
    {
      if (!SmootherPrototype::IsSetup())
        throw(Exceptions::RuntimeError("Setup has not been called"));

      Teuchos::ParameterList  ifpack2List;

      // Forward the InitialGuessIsZero option to Ifpack2
      //  TODO: It might be nice to switch back the internal
      //        "zero starting solution" option of the ifpack2 object prec_ to his
      //        initial value at the end but there is no way right now to get
      //        the current value of the "zero starting solution" in ifpack2.
      //        It's not really an issue, as prec_  can only be used by this method.
      if (ifpack2Type_ == "CHEBYSHEV") {
        ifpack2List.set("chebyshev: zero starting solution", InitialGuessIsZero);
      }
      else if (ifpack2Type_ == "RELAXATION") {
        ifpack2List.set("relaxation: zero starting solution", InitialGuessIsZero);
      }
      else if (ifpack2Type_ == "ILUT") {
        static int warning_only_once=0;
        if ((warning_only_once++) == 0)
          *out_ << "MueLu::Ifpack2Smoother::Apply(): Warning: ILUT as a smoother must solve correction equations but not implemented yet." << std::endl;
        // TODO: ILUT using correction equation should be implemented in ifpack2 directly
        //       I think that an option named "zero starting solution"
        //       is also appropriate for ILUT
      }
      else {
        // TODO: When https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c2 is done
        // we should remove the if/else/elseif and just test if this
        // option is supported by current ifpack2 preconditioner
        throw(Exceptions::RuntimeError("Ifpack2Smoother::Apply(): Ifpack2 preconditioner '"+ifpack2Type_+"' not supported"));
      }
      prec_->setParameters(ifpack2List);

      Tpetra::MultiVector<SC,LO,GO,NO> &tpX = Utils::MV2NonConstTpetraMV(X);
      Tpetra::MultiVector<SC,LO,GO,NO> const &tpB = Utils::MV2TpetraMV(B);

      prec_->apply(tpB,tpX);

    }

    //@}

    //! @name Utilities
    //@{

    void Print(std::string prefix) {
      Teuchos::OSTab tab(out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "Ifpack2Smoother::Print()" << std::endl;
      //prec_->Print(*out_);
    }

    RCP<SmootherPrototype> Copy()
    {
      return rcp(new Ifpack2Smoother(*this) );
    }

    void CopyParameters(RCP<SmootherPrototype> source)
    {
      RCP<Ifpack2Smoother> ifpack2Smoo = rcp_dynamic_cast<Ifpack2Smoother>(source);
      //TODO check if dynamic cast fails
      ifpack2Type_ = ifpack2Smoo->ifpack2Type_; //TODO: Get() methods
      prec_ = ifpack2Smoo->prec_;
      A_ = ifpack2Smoo->A_;
      overlap_ = ifpack2Smoo->overlap_;
      list_ = ifpack2Smoo->list_;
    }

    //@}

  }; //class Ifpack2Smoother

} //namespace MueLu

#define MUELU_IFPACK2_SMOOTHER_SHORT

#endif //ifdef HAVE_MUELU_IFPACK2

#endif //ifndef MUELU_IFPACK2_SMOOTHER_HPP
// Note: Ifpack2 may be able to accept directly MueLu matrix
