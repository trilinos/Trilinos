/*
 * MueLu_NoFactory.cpp
 *
 *  Created on: Sep 13, 2011
 *      Author: wiesner
 */

#include "MueLu_NoFactory.hpp"

namespace MueLu {

  RCP<const NoFactory> NoFactory::noFactory_ = Teuchos::null;

}
