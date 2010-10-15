#ifndef MUELU_TRANSPFACTORY_HPP
#define MUELU_TRANSPFACTORY_HPP

#include <iostream>
#include "MueLu_BaseFactory.hpp"

/*!
  @class TransPFactory class.
  @brief Factory for building Smoothed Aggregation prolongators.
*/

namespace MueLu {

class TransPFactory : public BaseFactory {
  inline friend std::ostream& operator<<(std::ostream& os, TransPFactory &factory);

  public:
    //@{ Constructors/Destructors.
    TransPFactory() {}

    virtual ~TransPFactory() {}
    //@}

}; //class TransPFactory

std::ostream& operator<<(std::ostream& os, TransPFactory &factory) {
  os << "Printing a TransPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#endif //ifndef MUELU_TRANSPFACTORY_HPP
