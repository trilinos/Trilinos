#ifndef MUELU_SAPFACTORY_HPP
#define MUELU_SAPFACTORY_HPP

#include <iostream>
#include "MueLu_BaseFactory.hpp"

/*!
  @class SaPFactory class.
  @brief Factory for building Smoothed Aggregation prolongators.
*/

namespace MueLu {

class SaPFactory : public BaseFactory {
  inline friend std::ostream& operator<<(std::ostream& os, SaPFactory &factory);

  public:
    //@{ Constructors/Destructors.
    SaPFactory() {}

    virtual ~SaPFactory() {}
    //@}

}; //class SaPFactory

std::ostream& operator<<(std::ostream& os, SaPFactory &factory) {
  os << "Printing an SaPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_SAPFACTORY_HPP
