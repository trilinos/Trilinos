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
    bool NoFactory::instanceFlag_ = false;
    NoFactory* NoFactory::UserDefined = NULL;

    NoFactory::ptrGetInstance UserDefinedVariable = NoFactory::get;
}

//MueLu::NoFactory::~NoFactory()


//static const MueLu::NoFactory* MueLu::NoFactory::get()

