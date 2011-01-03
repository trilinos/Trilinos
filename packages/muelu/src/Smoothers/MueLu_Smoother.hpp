#ifndef MUELU_SMOOTHER_HPP
#define MUELU_SMOOTHER_HPP

#include "MueLu_SmootherBase.hpp"
#include "MueLu_SmootherPrototype.hpp"

namespace MueLu {

template <class ScalarType,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class Level;

/*!
  @class Smoother
  @brief Smoother class. This is just a stub for a smoother.  It doesn't do any real computation.
*/

  template<class ScalarType,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class Smoother : public SmootherPrototype<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>
  {

#include "MueLu_UseShortNames.hpp"

  private:

    LO nIts_;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:

    //! Constructor
    Smoother() : out_(this->getOStream()) {std::cout << "Instantiating a new smoother" << std::endl;}

    virtual ~Smoother() {}

    void SetNIts(LO nIts) {
      nIts_ = nIts;
    }

    LO GetNIts() {
      return nIts_;
    }

    void CopyParameters(RCP<Smoother> source)
    {
      nIts_ = source.GetNIts();
    }

    void Setup(RCP<Level<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > level) {
      Teuchos::OSTab tab(out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "Smoother::Setup()" << std::endl;
      SmootherPrototype::IsSetup(true);
    }

    void Apply(RCP<MultiVector> x, RCP<MultiVector> const rhs, bool InitialGuessIsZero)
    {
      Teuchos::OSTab tab(out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "Smoother::Apply()" << std::endl;
    }

    void Print(std::string prefix) {
      Teuchos::OSTab tab(out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "Smoother::Print()" << std::endl;
    }

  }; //class Smoother

} //namespace MueLu

#define MUELU_SMOOTHER_SHORT

#endif //ifndef MUELU_SMOOTHER_HPP
