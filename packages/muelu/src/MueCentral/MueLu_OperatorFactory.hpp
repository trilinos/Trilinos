#ifndef MUELU_OPERATORFACTORY_HPP
#define MUELU_OPERATORFACTORY_HPP

#include <iostream>

#include "Teuchos_ParameterList.hpp"
#include "MueLu_Needs.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  /*!
    @class Basic class for operator factories (e.g., R, P, and A_coarse).
    @brief Basic class for factories that build operators.

    Derives from Needs, but with an additional virtual Build method.
  */

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class OperatorFactory : public Needs {

#include "MueLu_UseShortNames.hpp"
    
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
    virtual bool Build(Level & i, Level & j) const = 0;

    //@}

  }; //class OperatorFactory

} //namespace MueLu

#define MUELU_OPERATORFACTORY_SHORT

#endif //ifndef MUELU_OPERATORFACTORY_HPP
