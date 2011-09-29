/*
 * MueLu_NoFactory.cpp
 *
 *  Created on: Sep 13, 2011
 *      Author: wiesner
 */


#include "MueLu_NoFactory.hpp"

// static variables
namespace MueLu
{
RCP<NoFactory> NoFactory::UserDefined = Teuchos::null;

NoFactory::ptrGetInstance UserDefinedVariable = NoFactory::get;
}

//MueLu::NoFactory::~NoFactory()


//static const MueLu::NoFactory* MueLu::NoFactory::get()

