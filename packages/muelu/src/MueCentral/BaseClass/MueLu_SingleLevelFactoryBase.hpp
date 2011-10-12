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
    virtual void Build(Level & currentLevel) const = 0;

    //!
    virtual void CallBuild(Level & requestedLevel) const {
      Build(requestedLevel);
    }

    //!
    virtual void CallDeclareInput(Level & requestedLevel) const {
      DeclareInput(requestedLevel);
    }
 //@}

  }; //class SingleLevelFactoryBase

} //namespace MueLu

#define MUELU_SINGLELEVELFACTORY_SHORT
#endif //ifndef MUELU_SINGLELEVELFACTORY_HPP
