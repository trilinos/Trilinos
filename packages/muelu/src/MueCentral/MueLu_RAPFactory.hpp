#ifndef MUELU_RAPFACTORY_HPP
#define MUELU_RAPFACTORY_HPP

#include <iostream>
#include "MueLu_BaseFactory.hpp"

/*!
  @class RAPFactory class.
  @brief Factory for building coarse matrices.
*/

namespace MueLu {

class RAPFactory : public BaseFactory {

  inline friend std::ostream& operator<<(std::ostream& os, RAPFactory &factory);

  public:
    //@{ Constructors/Destructors.
    RAPFactory() {}

    virtual ~RAPFactory() {}
    //@}

    //@{ Set/Get methods.
    //@}


}; //class RAPFactory

std::ostream& operator<<(std::ostream& os, RAPFactory &factory) {
  os << "Printing RAPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_RAPFACTORY_HPP
