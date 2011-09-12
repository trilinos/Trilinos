#ifndef MUELU_IFPACK_SMOOTHER_HPP
#define MUELU_IFPACK_SMOOTHER_HPP

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_IFPACK

#include "Ifpack.h"

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  /*!
    @class IfpackSmoother
    @brief Class that encapsulates Ifpack smoothers.
    
    This class creates an Ifpack preconditioner factory.  The factory creates a smoother based on the
    type and ParameterList passed into the constructor.  See the constructor for more information.
  */
  class IfpackSmoother : public SmootherPrototype<double,int,int>
  {

    typedef double Scalar;
    typedef int    LocalOrdinal;
    typedef int    GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;
#include "MueLu_UseShortNames.hpp"

  public:

    //! @name Constructors / destructors
    //@{

    /*! @brief Constructor

        The options passed into IfpackSmoother are those given in the Ifpack user's manual.

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

        See also Ifpack_PointRelaxation, Ifpack_Chebyshev, Ifpack_ILU.
    */
    IfpackSmoother(std::string const & type, Teuchos::ParameterList const & paramList = Teuchos::ParameterList(), LO const &overlap=0) //TODO: empty paramList valid for Ifpack??
      : type_(type), paramList_(paramList), overlap_(overlap)
    { }

    //! Destructor
    virtual ~IfpackSmoother() {}

    //@}

    //! @name Set/Get methods

    //@{

    //! Set smoother parameters
    void SetParameters(Teuchos::ParameterList const & paramList) {
      paramList_ = paramList;

      if (SmootherPrototype::IsSetup()) {
        // It might be invalid to change parameters after the setup, but it depends entirely on Ifpack implementation.
        // TODO: I don't know if Ifpack returns an error code or exception or ignore parameters modification in this case...

        Teuchos::ParameterList nonConstParamList = paramList; // because Ifpack SetParameters() input argument is not const...
        prec_->SetParameters(nonConstParamList);
      }
    }

    //! Get smoother parameters
    Teuchos::ParameterList const & GetParameters() { return paramList_; }

    //JG: I'm not sure if it's a good idea to provide Get/Set NIts (for code maintainability)
    
    //     /*! @brief Set the number of smoothing sweeps/degree.
    //
    //        If the smoother is relaxation, this sets the number of sweeps.
    //        If the smoother is Chebyshev, this sets the polynomial degree.
    //
    //        Note:  This can be called after the preconditioner is set up, i.e., after
    //        calling IfpackSmoother::Setup().
    //     */
    //     void SetNIts(LO const &nIts) {
    //       if (!SmootherPrototype::IsSetup()) //FIXME precond doesn't have to be setup
    //         throw(Exceptions::RuntimeError("Call Setup before setting sweeps"));
    //       if (type_ == "point relaxation stand-alone") paramList_.set("relaxation: sweeps", nIts);
    //       else if (type_ == "Chebyshev")               paramList_.set("chebyshev: degree", nIts);
    //       else throw(Exceptions::RuntimeError("SetNIts: unknown smoother type"));
    //       prec_->SetParameters(paramList_);
    //     }
    //
    //     /*! @brief Get the number of smoothing sweeps.
    //
    //        If the smoother is relaxation, this returns the number of sweeps.
    //        If the smoother is Chebyshev, this returns the polynomial degree.
    //     */
    //     LO GetNIts() {
    //       if (type_ == "point relaxation stand-alone")
    //       {
    //         if (paramList_.isParameter("relaxation: sweeps") == false)
    //           throw(Exceptions::RuntimeError("number of iterations is not set"));
    //         return paramList_.get("relaxation: sweeps",1);
    //       } else if (type_ == "Chebyshev") {
    //         if (paramList_.isParameter("chebyshev: degree") == false)
    //           throw(Exceptions::RuntimeError("Chebyshev degree is not set"));
    //         return paramList_.get("chebyshev: degree",1);
    //       } else 
    //         throw(Exceptions::RuntimeError("GetNIts: unknown smoother type"));
    //     }

    //@}

    //! @name Computational methods.
    //@{

    /*! @brief Set up the smoother.

        This creates the underlying Ifpack smoother object, copies any parameter list options
        supplied to the constructor to the Ifpack object, and computes the preconditioner.
    */
    void Setup(Level &level) {
      TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == true, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Setup(): Setup() has already been called"); //TODO: Valid. To be replace by a warning.

      RCP<Teuchos::FancyOStream> out_ = this->getOStream();
      Teuchos::OSTab tab(out_);

      A_ = level.Get< RCP<Operator> >("A");

      // output information
      std::ostringstream buf; buf << level.GetLevelID();
      std::string prefix = "Smoother (level " + buf.str() + ") : ";
      LO rootRank = out_->getOutputToRootOnly();
      out_->setOutputToRootOnly(0);
      *out_ << prefix << "# global rows = " << A_->getGlobalNumRows()
            << ", estim. global nnz = " << A_->getGlobalNumEntries() << std::endl;

      RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);
      Ifpack factory;
      prec_ = rcp(factory.Create(type_, &(*epA), overlap_));
      if (type_ == "Chebyshev") {
        Scalar maxEigenValue = paramList_.get("chebyshev: max eigenvalue",(Scalar)-1.0);
        if (maxEigenValue == -1.0) {
          maxEigenValue = Utils::PowerMethod(*A_,true,10,1e-4);
          paramList_.set("chebyshev: max eigenvalue",maxEigenValue);
        }
        *out_ << prefix << "Ifpack Chebyshev, degree " << paramList_.get("chebyshev: degree",1) << std::endl;
        *out_ << prefix << "lambda_min=" << paramList_.get("chebyshev: min eigenvalue",-1.0)
              << ", lambda_max=" << paramList_.get("chebyshev: max eigenvalue",-1.0) << std::endl;
      }
      out_->setOutputToRootOnly(rootRank);
      prec_->SetParameters(paramList_);
      prec_->Compute();

      SmootherPrototype::IsSetup(true);
    }

    /*! @brief Apply the preconditioner.

        Solves the linear system <tt>AX=B</tt> using the constructed smoother.

        @param X initial guess
        @param B right-hand side
        @param InitialGuessIsZero (optional) If false, some work can be avoided.  Whether this actually saves any work depends on the underlying Ifpack implementation.
    */
    void Apply(MultiVector& X, MultiVector const &B, bool const &InitialGuessIsZero=false) const {
      TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Apply(): Setup() has not been called");

      RCP<Teuchos::FancyOStream> out_ = this->getOStream();

      // Forward the InitialGuessIsZero option to Ifpack
      Teuchos::ParameterList  paramList;
      if (type_ == "Chebyshev") {
        paramList.set("chebyshev: zero starting solution", InitialGuessIsZero);
      } else if (type_ == "point relaxation stand-alone") {
        paramList.set("relaxation: zero starting solution", InitialGuessIsZero);
      } else if  (type_ == "ILU") {
         if (InitialGuessIsZero == false) {
           Teuchos::OSTab tab(out_);
           *out_ << "WARNING:  Ifpack2 ILUT has no provision for a nonzero initial guess." << std::endl;
         }
      } else {
        // TODO: When https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c2 is done
        // we should remove the if/else/elseif and just test if this
        // option is supported by current ifpack2 preconditioner
        TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,"IfpackSmoother::Apply(): Ifpack preconditioner '"+type_+"' not supported");
      }
      prec_->SetParameters(paramList);

      // Apply
      Epetra_MultiVector &epX = Utils::MV2NonConstEpetraMV(X);
      Epetra_MultiVector const &epB = Utils::MV2EpetraMV(B);
      prec_->ApplyInverse(epB, epX);
    }

    //@}

    //! @name Utilities
    //@{

    RCP<SmootherPrototype> Copy() const {
      return rcp(new IfpackSmoother(*this) );
    }

    //@}

  private:

    //! ifpack-specific key phrase that denote smoother type
    std::string type_;

    //! parameter list that is used by Ifpack internally
    Teuchos::ParameterList paramList_;

    //! overlap when using the smoother in additive Schwarz mode
    LO overlap_;

    //! pointer to Ifpack solver object
    RCP<Ifpack_Preconditioner> prec_;

    //! Operator. Not used directly, but held inside of prec_. So we have to keep an RCP pointer to it!
    RCP<Operator> A_;

  }; // class IfpackSmoother

} // namespace MueLu

#define MUELU_IFPACK_SMOOTHER_SHORT
#endif //ifdef HAVE_MUELU_IFPACK
#endif //ifndef MUELU_IFPACK_SMOOTHER_HPP

// For describe()
//       std::string label;
//       if (type == "point relaxation stand-alone")
//         label = "Ifpack: " + paramList_.get("relaxation: type","unknown relaxation");
//       else
//         label = "Ifpack: " + type;
//       //TODO      SmootherBase::SetType(label);
