/*
 * MueLu_NoFactory.hpp
 *
 *  Created on: Sep 13, 2011
 *      Author: wiesner
 */

#ifndef MUELU_NOFACTORY_HPP_
#define MUELU_NOFACTORY_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu
{
class Level;

/*!
      @class NoFactory class.
      @brief NoFactory that is used for data stored in level class for that no generating factory is available/necessary.
 */
class NoFactory : public SingleLevelFactoryBase {
private:
  static RCP<NoFactory> UserDefined; // static NoFactory instance for user defined "factories"

  //! Constructor.
  NoFactory() {}

public:
  //@{ Destructor.


  //! Destructor.
  virtual ~NoFactory() { /*NoFactory::UserDefined = Teuchos::null;*/ }
  //@}

  //! Input
  //@{

  void DeclareInput(Level &currentLevel) const { }

  //@}

  //@{
  //! @name Build methods.

  //! Build. The NoFactory has no support for a Build routine. throw error
  void Build(Level & currentLevel) const
  {
    TEST_FOR_EXCEPTION(1, MueLu::Exceptions::RuntimeError, "NoFactory::Build(): the Build method of NoFactory cannot be called.");
  }

  //!
  void NewBuild(Level & requestedLevel) const {
    return Build(requestedLevel);
  }

  //!
  void callDeclareInput(Level & requestedLevel) const {
    DeclareInput(requestedLevel); }

  static const MueLu::NoFactory* get() {
    if(NoFactory::UserDefined == Teuchos::null)
    {
      NoFactory::UserDefined = rcp(new NoFactory());
      return NoFactory::UserDefined.get();
    }
    else
      return MueLu::NoFactory::UserDefined.get();
  }

  typedef const NoFactory* (*ptrGetInstance)();

}; // end NoFactory

} // end namespace MueLu

#endif
