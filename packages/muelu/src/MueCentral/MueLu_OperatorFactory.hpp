#ifndef MUELU_OPERATORFACTORY_HPP
#define MUELU_OPERATORFACTORY_HPP

#include <iostream>

#include "Teuchos_ParameterList.hpp"
#include "MueLu_BaseFactory.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

/*!
  @class Basic class for operator factories (e.g., R, P, and A_coarse).
  @brief Basic class for factories that build operators.

  Very similar to BaseFactory, but with an additional virtual Build method.
*/

template<class Scalar, class LO, class GO, class Node>
class OperatorFactory : public BaseFactory {

  private:

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    OperatorFactory() {}

    //! Destructor.
    virtual ~OperatorFactory() {}
    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    virtual bool Build(Level<Scalar, LO, GO, Node> &i,Level<Scalar, LO, GO, Node> &j) = 0;

    //@}

}; //class OperatorFactory

} //namespace MueLu

#endif //ifndef MUELU_OPERATORFACTORY_HPP
