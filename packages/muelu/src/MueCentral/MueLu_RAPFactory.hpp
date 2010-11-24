#ifndef MUELU_RAPFACTORY_HPP
#define MUELU_RAPFACTORY_HPP

#include <iostream>
#include "MueLu_OperatorFactory.hpp"

namespace MueLu {
/*!
  @class RAPFactory class.
  @brief Factory for building coarse matrices.
*/
  template<class Scalar, class LO, class GO, class NO, class LMO>
class RAPFactory : public OperatorFactory<Scalar,LO,GO,NO,LMO> {

  typedef Level<Scalar,LO,GO,NO,LMO> Level;

  //JG to JJH: use Teuchos::Describable instead ?
  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, RAPFactory<AA,BB,CC,DD,EE> &factory);

  public:
    //@{ Constructors/Destructors.
    RAPFactory() {}

    virtual ~RAPFactory() {}
    //@}

    //@{ Build methods.
    bool Build(Level &fineLevel, Level &coarseLevel) {
     Teuchos::OSTab tab(this->out_);  *(this->out_) << "RAPFactory: Building a coarse operator" << std::endl; return true;}
    //@}


}; //class RAPFactory

//! Friend print method.
  template<class Scalar, class LO, class GO, class NO, class LMO>
std::ostream& operator<<(std::ostream& os, RAPFactory<Scalar,LO,GO,NO,LMO> &factory) {
  os << "Printing RAPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_RAPFACTORY_HPP
