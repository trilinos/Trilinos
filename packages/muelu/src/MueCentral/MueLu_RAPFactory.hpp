#ifndef MUELU_RAPFACTORY_HPP
#define MUELU_RAPFACTORY_HPP

#include <iostream>
#include "MueLu_OperatorFactory.hpp"


namespace MueLu {
/*!
  @class RAPFactory class.
  @brief Factory for building coarse matrices.
*/
template<class Scalar, class LO, class GO, class Node>
class RAPFactory : public OperatorFactory<Scalar,LO,GO,Node> {

  template<class AA, class BB, class CC, class DD>
  inline friend std::ostream& operator<<(std::ostream& os, RAPFactory<AA,BB,CC,DD> &factory);

  public:
    //@{ Constructors/Destructors.
    RAPFactory() {}

    virtual ~RAPFactory() {}
    //@}

    //@{ Build methods.
    bool Build(Level<Scalar,LO,GO,Node> &fineLevel, Level<Scalar,LO,GO,Node> &coarseLevel) {
     Teuchos::OSTab tab(this->out_);  *(this->out_) << "RAPFactory: Building a coarse operator" << std::endl; return true;}
    //@}


}; //class RAPFactory

//! Friend print method.
template<class Scalar, class LO, class GO, class Node>
std::ostream& operator<<(std::ostream& os, RAPFactory<Scalar,LO,GO,Node> &factory) {
  os << "Printing RAPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_RAPFACTORY_HPP
