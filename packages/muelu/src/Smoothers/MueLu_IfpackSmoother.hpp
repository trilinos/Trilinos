#ifndef MUELU_IFPACK_SMOOTHER_HPP
#define MUELU_IFPACK_SMOOTHER_HPP

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Utilities.hpp"

//#ifdef HAVE_MUELU_IFPACK //FIXME doesn't work right now
#include "Ifpack.h"

namespace MueLu {

template <class ScalarType,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class Level;

/*!
  @class IfpackSmoother
  @brief Class that encapsulates Ifpack smoothers.

  It does this by invoking the Ifpack preconditioner factory.
*/

  template<class ScalarType,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class IfpackSmoother : public SmootherPrototype<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

#include "MueLu_UseShortNames.hpp"

  private:

    //! Ifpack-specified key phrase that denote smoother type
    std::string type_;
    LO nIts_;
    LO overlap_;
    //RCP<Ifpack_Preconditioner> prec_;
    Ifpack_Preconditioner* prec_;
    Teuchos::RCP<Operator> A_;
    Teuchos::ParameterList list_;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:

    //! @name Constructors / destructors
    //@{
    IfpackSmoother(std::string const & type, Teuchos::ParameterList & list)
      : type_(type), list_(list), out_(this->getOStream())
    {
      MueLu_cout(Teuchos::VERB_HIGH) << "Instantiating a new smoother" << std::endl;
      overlap_ = list.get("Overlap",(LO) 0);
    }

    virtual ~IfpackSmoother() {}
    //@}

    void SetNIts(LO nIts) {
      nIts_ = nIts;
    }

    LO GetNIts() {
      return nIts_;
    }

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

//#endif //ifdef HAVE_MUELU_IFPACK //FIXME doesn't work right now

#endif //ifndef MUELU_IFPACK_SMOOTHER_HPP
