#ifndef MUELU_IFPACK_SMOOTHER_HPP
#define MUELU_IFPACK_SMOOTHER_HPP

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_IFPACK

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities.hpp"

#include "Ifpack.h"

namespace MueLu {

class Level;

/*!
  @class IfpackSmoother
  @brief Class that encapsulates Ifpack smoothers.

  This class creates an Ifpack preconditioner factory.  The factory creates a smoother based on the
  type and ParameterList passed into the constructor.  See the constructor for more information.
*/

  template<class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class IfpackSmoother : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

#include "MueLu_UseShortNames.hpp"

  private:

    //! Ifpack-specific key phrase that denote smoother type (not to be confused with SmootherBase::Type_)
    std::string ifpackType_;
    //! overlap when using the smoother in additive Schwarz mode
    LO overlap_;
    RCP<Ifpack_Preconditioner> prec_;
    //! matrix operator 
    RCP<Operator> A_;
    //! parameter list that is used by Ifpack internally
    Teuchos::ParameterList list_;

  protected:
    RCP<Teuchos::FancyOStream> out_;

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
    IfpackSmoother(std::string const & type, Teuchos::ParameterList & list)
      : ifpackType_(type), list_(list), out_(this->getOStream())
    {
      overlap_ = list.get("Overlap",(LO) 0);
      std::string label;
      if (type == "point relaxation stand-alone")
        label = "Ifpack: " + list.get("relaxation: type","unknown relaxation");
      else
        label = "Ifpack: " + type;
      SmootherBase::SetType(label);
      SmootherPrototype::IsSetup(false);
    }

    //! Destructor
    virtual ~IfpackSmoother() {}
    //@}

    //! @name Set/Get methods

    //@{

    /*! @brief Set the number of smoothing sweeps/degree.

       If the smoother is relaxation, this sets the number of sweeps.
       If the smoother is Chebyshev, this sets the polynomial degree.

       Note:  This can be called after the preconditioner is set up, i.e., after
       calling IfpackSmoother::Setup().
    */
    void SetNIts(LO const &nIts) {
      if (!SmootherPrototype::IsSetup()) //FIXME precond doesn't have to be setup
        throw(Exceptions::RuntimeError("Call Setup before setting sweeps"));
      if (ifpackType_ == "point relaxation stand-alone") list_.set("relaxation: sweeps", nIts);
      else if (ifpackType_ == "Chebyshev")               list_.set("chebyshev: degree", nIts);
      else throw(Exceptions::RuntimeError("SetNIts: unknown smoother type"));
      prec_->SetParameters(list_);
    }

    /*! @brief Get the number of smoothing sweeps.

       If the smoother is relaxation, this returns the number of sweeps.
       If the smoother is Chebyshev, this returns the polynomial degree.
    */
    LO GetNIts() {
      if (ifpackType_ == "point relaxation stand-alone")
      {
        if (list_.isParameter("relaxation: sweeps") == false)
          throw(Exceptions::RuntimeError("number of iterations is not set"));
        return list_.get("relaxation: sweeps",1);
      } else if (ifpackType_ == "Chebyshev") {
        if (list_.isParameter("chebyshev: degree") == false)
          throw(Exceptions::RuntimeError("Chebyshev degree is not set"));
        return list_.get("chebyshev: degree",1);
      } else 
        throw(Exceptions::RuntimeError("GetNIts: unknown smoother type"));
    }
    //@}

    //! @name Computational methods.
    //@{

    /*! @brief Set up the smoother.

        This creates the underlying Ifpack smoother object, copies any parameter list options
        supplied to the constructor to the Ifpack object, and computes the preconditioner.
    */
    void Setup(Level &level) {
      Teuchos::OSTab tab(out_);

      A_ = level.Get< RCP<Operator> >("A");
      RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);
      Ifpack factory;
      prec_ = rcp(factory.Create(ifpackType_, &(*epA), overlap_));
      prec_->SetParameters(list_);
      prec_->Compute();

      SmootherPrototype::IsSetup(true);
    }


    /*! @brief Apply the preconditioner.

        Solves the linear system <tt>AX=B</tt> using the constructed smoother.

        @param X initial guess
        @param B right-hand side
        @param InitialGuessIsZero (optional) If false, some work can be avoided.  Whether this actually saves any work depends on the underlying Ifpack implementation.
    */
    void Apply(MultiVector& X, MultiVector const &B, bool const &InitialGuessIsZero=false)
    {
      if (!SmootherPrototype::IsSetup())
        throw(Exceptions::RuntimeError("Setup has not been called"));
      Teuchos::ParameterList  ifpackList;
      ifpackList.set("relaxation: zero starting solution", InitialGuessIsZero);
      prec_->SetParameters(ifpackList);

      Epetra_MultiVector &epX = Utils::MV2NonConstEpetraMV(X);
      Epetra_MultiVector const &epB = Utils::MV2EpetraMV(B);

      prec_->ApplyInverse(epB,epX);
    }

    //@}

    //! @name Utilities
    //@{

    void Print(std::string prefix) const {
      Teuchos::OSTab tab(out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "IfpackSmoother::Print()" << std::endl;
      prec_->Print(*out_);
    }

    RCP<SmootherPrototype> Copy() const
    {
      return rcp(new IfpackSmoother(*this) );
    }

    void CopyParameters(RCP<SmootherPrototype> source)
    {
      RCP<IfpackSmoother> ifpackSmoo = rcp_dynamic_cast<IfpackSmoother>(source);
      //TODO check if dynamic cast fails
      ifpackType_ = ifpackSmoo->ifpackType_;
      prec_ = ifpackSmoo->prec_;
      A_ = ifpackSmoo->A_;
      overlap_ = ifpackSmoo->overlap_;
      list_ = ifpackSmoo->list_;
    }

    //@}

  }; //class IfpackSmoother

} //namespace MueLu

#define MUELU_IFPACK_SMOOTHER_SHORT

#endif //ifdef HAVE_MUELU_IFPACK

#endif //ifndef MUELU_IFPACK_SMOOTHER_HPP
