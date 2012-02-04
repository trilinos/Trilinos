#ifndef MUELU_TWOLEVELFACTORY_HPP
#define MUELU_TWOLEVELFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  /*!
    @class TwoLevelFactoryBase class.
    @brief Base class for factories that use two levels (fineLevel and coarseLevel).

    Examples of such factories are R, P, and A_coarse.

    @ingroup MueLuBaseClasses
  */


  class TwoLevelFactoryBase : public FactoryBase {

  public:

    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    TwoLevelFactoryBase() {}

    //! Destructor.
    virtual ~TwoLevelFactoryBase() {}

    //@}

    //! Input
    //@{

    /*! @brief Specifies the data that this class needs, and the factories that generate that data.

        If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
        will fall back to the settings in FactoryManager.
    */
    virtual void DeclareInput(Level &fineLevel, Level &coarseLevel) const = 0;

    //!
    virtual void CallDeclareInput(Level & requestedLevel) const {
      TEUCHOS_TEST_FOR_EXCEPTION(requestedLevel.GetPreviousLevel() == Teuchos::null, Exceptions::RuntimeError, "LevelID = " << requestedLevel.GetLevelID());
      DeclareInput(*requestedLevel.GetPreviousLevel(), requestedLevel);
    }

    //@}

    //! @name Build methods.
    //@{


    //! Build an object with this factory.
    virtual void Build(Level & fineLevel, Level & coarseLevel) const = 0;

    //!
    virtual void CallBuild(Level & requestedLevel) const {
      TEUCHOS_TEST_FOR_EXCEPTION(requestedLevel.GetPreviousLevel() == Teuchos::null, Exceptions::RuntimeError, "LevelID = " << requestedLevel.GetLevelID());
      Build(*requestedLevel.GetPreviousLevel(), requestedLevel);
    }

    //@}

  }; //class TwoLevelFactoryBase


} //namespace MueLu

#define MUELU_TWOLEVELFACTORY_SHORT
#endif //ifndef MUELU_TWOLEVELFACTORY_HPP
