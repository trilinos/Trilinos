#ifndef MUELU_SAPFACTORY_HPP
#define MUELU_SAPFACTORY_HPP

#include <iostream>
#include "MueLu_OperatorFactory.hpp"

/*!
  @class SaPFactory class.
  @brief Factory for building Smoothed Aggregation prolongators.
*/

namespace MueLu {

template<class Scalar, class LO, class GO, class Node>
class SaPFactory : public OperatorFactory<Scalar,LO,GO,Node> {
  template<class AA, class BB, class CC, class DD>
  inline friend std::ostream& operator<<(std::ostream& os, SaPFactory<AA,BB,CC,DD> &factory);

  public:
    //@{ Constructors/Destructors.
    SaPFactory() {}

    virtual ~SaPFactory() {}
    //@}

    //@{ Build methods.
    bool Build(Level<Scalar,LO,GO,Node> &fineLevel, Level<Scalar,LO,GO,Node> &coarseLevel) {std::cout << "SaPFactory: Building a prolongator" << std::endl; return true;}
    //@}

}; //class SaPFactory

template<class Scalar, class LO, class GO, class Node>
std::ostream& operator<<(std::ostream& os, SaPFactory<Scalar,LO,GO,Node> &factory) {
  os << "Printing an SaPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_SAPFACTORY_HPP
