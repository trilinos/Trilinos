#ifndef MUELU_TRANSPFACTORY_HPP
#define MUELU_TRANSPFACTORY_HPP

#include <iostream>
#include "MueLu_OperatorFactory.hpp"

/*!
  @class TransPFactory class.
  @brief Factory for building Smoothed Aggregation prolongators.
*/

namespace MueLu {

template<class Scalar, class LO, class GO, class Node>
class TransPFactory : public OperatorFactory<Scalar,LO,GO,Node> {
  template<class AA, class BB, class CC, class DD>
  inline friend std::ostream& operator<<(std::ostream& os, TransPFactory<AA,BB,CC,DD> &factory);

  public:
    //@{ Constructors/Destructors.
    TransPFactory() {}

    virtual ~TransPFactory() {}
    //@}

    //@{ Build methods.
    bool Build(Level<Scalar,LO,GO,Node> &fineLevel, Level<Scalar,LO,GO,Node> &coarseLevel) {*(this->out_) << "TransPFactory: Building a restriction operator" << std::endl; return true;}
    //@}

}; //class TransPFactory

template<class Scalar, class LO, class GO, class Node>
std::ostream& operator<<(std::ostream& os, TransPFactory<Scalar,LO,GO,Node> &factory) {
  os << "Printing a TransPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_TRANSPFACTORY_HPP
