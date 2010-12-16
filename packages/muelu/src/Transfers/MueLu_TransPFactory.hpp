#ifndef MUELU_TRANSPFACTORY_HPP
#define MUELU_TRANSPFACTORY_HPP

#include <iostream>
#include "MueLu_OperatorFactory.hpp"

namespace MueLu {

  /*!
    @class TransPFactory class.
    @brief Factory for building restriction operators.
  */

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class TransPFactory : public OperatorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

    template<class AA, class BB, class CC, class DD, class EE>
    inline friend std::ostream& operator<<(std::ostream& os, TransPFactory<AA,BB,CC,DD,EE> &factory);

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    TransPFactory() {
      Teuchos::OSTab tab(this->out_);
      MueLu_cout(Teuchos::VERB_HIGH) << "TransPFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~TransPFactory() {}
    //@}

    //@{ Build methods.
    bool Build(Level & fineLevel, Level & coarseLevel) {
      Teuchos::OSTab tab(this->out_); MueLu_cout(Teuchos::VERB_HIGH) << "TransPFactory: Building a restriction operator" << std::endl; return true;
    }
    //@}

  }; //class TransPFactory

  //! Friend print function.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::ostream& operator<<(std::ostream& os, TransPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &factory) {
    os << "Printing a TransPFactory object" << std::endl;
    return os;
  }

} //namespace MueLu

#define MUELU_TRANSPFACTORY_SHORT

#endif //ifndef MUELU_TRANSPFACTORY_HPP
