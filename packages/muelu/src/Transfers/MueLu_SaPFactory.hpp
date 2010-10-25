#ifndef MUELU_SAPFACTORY_HPP
#define MUELU_SAPFACTORY_HPP

#include <iostream>
#include "MueLu_OperatorFactory.hpp"

namespace MueLu {

/*!
  @class SaPFactory class.
  @brief Factory for building Smoothed Aggregation prolongators.
*/

template<class Scalar, class LO, class GO, class Node>
class SaPFactory : public OperatorFactory<Scalar,LO,GO,Node> {
  template<class AA, class BB, class CC, class DD>
  inline friend std::ostream& operator<<(std::ostream& os, SaPFactory<AA,BB,CC,DD> &factory);

  public:
    //@{ Constructors/Destructors.
    //! Constructor.
    SaPFactory() {}

    //! Destructor.
    virtual ~SaPFactory() {}
    //@}

    //@{ Build methods.

    //! Build method.
    bool Build(Level<Scalar,LO,GO,Node> &fineLevel, Level<Scalar,LO,GO,Node> &coarseLevel) {
      Teuchos::OSTab tab(this->out_); *(this->out_) << "SaPFactory: Building a prolongator" << std::endl; return true;
    }
    //@}

}; //class SaPFactory

//! Friend print function.
template<class Scalar, class LO, class GO, class Node>
std::ostream& operator<<(std::ostream& os, SaPFactory<Scalar,LO,GO,Node> &factory) {
  os << "Printing an SaPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_SAPFACTORY_HPP
