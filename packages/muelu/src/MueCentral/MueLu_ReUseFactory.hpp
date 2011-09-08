#ifndef MUELU_REUSEFACTORY_HPP
#define MUELU_REUSEFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_FactoryBase.hpp"

namespace MueLu {

  class ReUseFactory : public FactoryBase {

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ReUseFactory() //TODO Add argument (std::string, factory) to be able to check if data available
    {

    }

    //! Destructor.
    virtual ~ReUseFactory() {}
    //@}

    //! @name Set/Get methods.
    //@{

    //TODO Add option to do the build if not available or to throw exception

    //@}

    //! Input
    //@{

    // TAW: no DeclareInput, since it derives from FactoryBase???
    //void DeclareInput(Level &currentLevel) const {
      // TODO
    //}

    //@}

    //! @name Build methods.
    //@{

    /*! @brief Build aggregates. */
    virtual bool NewBuild(Level & requestedLevel) const {
      // TODO: add some tests here, ie throw exception if reuse not possible and option ReBuild=false
      // or, better, call currentLevel->Get(string_,factory_); to transfer data to "this"

      return true; // ??
    }

    virtual void callDeclareInput(Level & requestedLevel) const {
      return;
    }

  }; // class ReUseFactory

} //namespace MueLu

#define MUELU_REUSEFACTORY_SHORT
#endif //ifndef MUELU_REUSEFACTORY_HPP
