#ifndef MUELU_SINGLELEVELFACTORY_HPP
#define MUELU_SINGLELEVELFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryBase.hpp"

namespace MueLu {
  class Level;

  /*!
    @class Base class for factories (e.g., Aggregation, Smoothers).
    @brief Base class for factories that use one levels (currentLevel).
  */
  class SingleLevelFactoryBase : public FactoryBase {

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    SingleLevelFactoryBase() {}

    //! Destructor.
    virtual ~SingleLevelFactoryBase() {}
    //@}

    //! Input
    //@{

    virtual void DeclareInput(Level &currentLevel) const = 0;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    virtual bool Build(Level & currentLevel) const = 0;

    //!
    virtual bool NewBuild(Level & requestedLevel) const {
      return Build(requestedLevel);
    }
 //@}

  }; //class SingleLevelFactoryBase

} //namespace MueLu

#define MUELU_SINGLELEVELFACTORY_SHORT
#endif //ifndef MUELU_SINGLELEVELFACTORY_HPP
