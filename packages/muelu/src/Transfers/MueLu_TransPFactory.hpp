#ifndef MUELU_TRANSPFACTORY_HPP
#define MUELU_TRANSPFACTORY_HPP

#include <iostream>
#include "MueLu_OperatorFactory.hpp"

namespace MueLu {

  /*!
    @class TransPFactory class.
    @brief Factory for building restriction operators.
  */

  template<class Scalar, class LO, class GO, class NO, class LMO>
  class TransPFactory : public OperatorFactory<Scalar,LO,GO,NO, LMO> {
    template<class AA, class BB, class CC, class DD, class EE>
    inline friend std::ostream& operator<<(std::ostream& os, TransPFactory<AA,BB,CC,DD,EE> &factory);

    typedef Level<Scalar,LO,GO,NO, LMO> Level;

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
  template<class Scalar, class LO, class GO, class NO, class LMO>
  std::ostream& operator<<(std::ostream& os, TransPFactory<Scalar,LO,GO,NO, LMO> &factory) {
    os << "Printing a TransPFactory object" << std::endl;
    return os;
  }

} //namespace MueLu

#endif //ifndef MUELU_TRANSPFACTORY_HPP
