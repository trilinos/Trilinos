#ifndef MUELU_IFPACK_SMOOTHER_HPP
#define MUELU_IFPACK_SMOOTHER_HPP

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities.hpp"

#ifdef HAVE_MUELU_IFPACK
#include "Ifpack.h"

namespace MueLu {

template <class ScalarType,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class Level;

/*!
  @class IfpackSmoother
  @brief Class that encapsulates Ifpack smoothers.

  This class creates an Ifpack preconditioner factory.  The factory creates a smoother based on the
  type and ParameterList passed into the constructor.  See the constructor for more information.
*/

  template<class ScalarType,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class IfpackSmoother : public SmootherPrototype<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

#include "MueLu_UseShortNames.hpp"

  private:

    //! Ifpack-specified key phrase that denote smoother type
    std::string type_;
    //! number of smoother sweeps
    LO nIts_;
    //! overlap when using the smoother in additive Schwarz mode
    LO overlap_;
    //RCP<Ifpack_Preconditioner> prec_;
    Ifpack_Preconditioner* prec_;
    //! matrix operator 
    Teuchos::RCP<Operator> A_;
    //! parameter list that is used by Ifpack internally
    Teuchos::ParameterList list_;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

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
         - ILU
            - <tt>type</tt> = <tt>ILU</tt>
            - parameter list options
                - <tt>fact: level-of-fill</tt>
    */
    IfpackSmoother(std::string const & type, Teuchos::ParameterList & list)
      : type_(type), list_(list), out_(this->getOStream())
    {
      MueLu_cout(Teuchos::VERB_HIGH) << "Instantiating a new smoother" << std::endl;
      overlap_ = list.get("Overlap",(LO) 0);
    }

    //! Destructor
    virtual ~IfpackSmoother() {}
    //@}

    //! @name Set/Get methods

    //@{
    void SetNIts(LO nIts) {
      nIts_ = nIts;
    }

    LO GetNIts() {
      return nIts_;
    }
    //@}

    void Setup(RCP<Level> level) {
      Teuchos::OSTab tab(out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "IfpackSmoother::Setup()" << std::endl;
      SmootherPrototype::IsSetup(true);
      A_ = level->GetA();
      RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);
      Ifpack factory;
      prec_ = factory.Create(type_, &(*epA), overlap_);
      prec_->SetParameters(list_);
      prec_->Compute();
    }

    void Apply(RCP<MultiVector> x, RCP<MultiVector> const rhs, bool InitialGuessIsZero)
    {
      if (InitialGuessIsZero)
        throw(Exceptions::NotImplemented("No logic for handling zero initial guesses"));
      Teuchos::OSTab tab(out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "IfpackSmoother::Apply()" << std::endl;
      RCP<Epetra_MultiVector> epX = Utils::MV2NonConstEpetraMV(x);
      RCP<const Epetra_MultiVector> epRhs = Utils::MV2NonConstEpetraMV(rhs);
      prec_->ApplyInverse(*epRhs,*epX);
    }

    void Print(std::string prefix) {
      Teuchos::OSTab tab(out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "IfpackSmoother::Print()" << std::endl;
      prec_->Print(*out_);
    }

    RCP<SmootherPrototype> Copy()
    {
      return rcp(new IfpackSmoother(*this) );
    }

    void CopyParameters(RCP<SmootherPrototype> source)
    {
      RCP<IfpackSmoother> ifpackSmoo = rcp_dynamic_cast<IfpackSmoother>(source);
      type_ = ifpackSmoo->type_;
      nIts_ = ifpackSmoo->nIts_;
      prec_ = ifpackSmoo->prec_;
      A_ = ifpackSmoo->A_;
    }

  }; //class IfpackSmoother

} //namespace MueLu

#define MUELU_IFPACK_SMOOTHER_SHORT

#endif //ifdef HAVE_MUELU_IFPACK

#endif //ifndef MUELU_IFPACK_SMOOTHER_HPP
