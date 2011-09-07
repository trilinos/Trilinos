#ifndef MUELU_SMOOTHERFACTORYBASE_HPP
#define MUELU_SMOOTHERFACTORYBASE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Types.hpp"

namespace MueLu {
  class Level;

  /*!
    @class Smoother factory base class.
    @brief Base class for smoother factories.

    This interface inherits from SingleLevelFactoryBase and defines an additional Build method (BuildSmoother)
    This new build method is used by Hierarchy for the CoarsestLevel.
  */

  class SmootherFactoryBase : public SingleLevelFactoryBase {

  public:
    //@{ Constructors/Destructors.
    SmootherFactoryBase() {}

    virtual ~SmootherFactoryBase() {}
    //@}

    //! @name Build methods.
    //@{

    //! Build pre-smoother and/or post-smoother
    virtual bool Build(Level & currentLevel) const = 0;
  
    virtual bool BuildSmoother(Level & currentLevel, PreOrPost const preOrPost = BOTH) const = 0;
    //@}

  }; //class SmootherFactoryBase

} //namespace MueLu

#define MUELU_SMOOTHERFACTORYBASE_SHORT

#endif //ifndef MUELU_SMOOTHERFACTORYBASE_HPP

//TODO: remove this interface?
