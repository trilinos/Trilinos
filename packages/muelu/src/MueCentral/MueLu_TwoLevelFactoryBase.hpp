#ifndef MUELU_TWOLEVELFACTORY_HPP
#define MUELU_TWOLEVELFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {
  class Level;

  /*!
    @class Base class for factories (e.g., R, P, and A_coarse).
    @brief Base class for factories that use two levels (fineLevel and coarseLevel).
  */

  class TwoLevelFactoryBase : public FactoryBase {

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    TwoLevelFactoryBase() {}

    //! Destructor.
    virtual ~TwoLevelFactoryBase() {}
    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    virtual bool Build(Level & fineLevel, Level & coarseLevel) const = 0;

    //!
    virtual bool NewBuild(Level & requestedLevel) const {
      return Build(*requestedLevel.GetPreviousLevel(), requestedLevel);
    }

    //@}

  }; //class TwoLevelFactoryBase

} //namespace MueLu

#define MUELU_TWOLEVELFACTORY_SHORT
#endif //ifndef MUELU_TWOLEVELFACTORY_HPP
