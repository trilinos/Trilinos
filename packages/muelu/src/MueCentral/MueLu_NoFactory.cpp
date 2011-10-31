/*
 * MueLu_NoFactory.cpp
 *
 *  Created on: Sep 13, 2011
 *      Author: wiesner
 */

#include "MueLu_NoFactory.hpp"

namespace MueLu {

  NoFactory::NoFactory() { }

  NoFactory::~NoFactory() { }

  void NoFactory::CallBuild(Level & requestedLevel) const {  
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::NoFactory::Build(): this method cannot be called.");
  }

  void NoFactory::CallDeclareInput(Level & requestedLevel) const {  }

  const RCP<const NoFactory> NoFactory::getRCP() {
    if(noFactory_ == Teuchos::null) {
      noFactory_ = rcp(new NoFactory());
    }

    return noFactory_;
  }

  const NoFactory* NoFactory::get() {
    return getRCP().get();
  }

  RCP<const NoFactory> NoFactory::noFactory_ = Teuchos::null;

} // namespace MueLu
