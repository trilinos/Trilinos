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

  template<class Scalar, class LO, class GO, class NO, class LMO>
  class OperatorFactory : public Needs {
    
    typedef Level<Scalar, LO, GO, NO, LMO> Level;

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
    virtual bool Build(Level & i, Level & j) = 0;

    //@}

  }; //class OperatorFactory

} //namespace MueLu

#endif //ifndef MUELU_OPERATORFACTORY_HPP
