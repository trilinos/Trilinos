// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EVIL_BASE_DEF_HPP
#define EVIL_BASE_DEF_HPP


#include "EvilBase_decl.hpp"

// Include the subclasses that we are going to instantiate in the factory
// (i.e. evil)!  NOTE: We need to include the possible function definitions
// here in case we are doing implicit instantiation!
#include "AEvil.hpp"
#include "BEvil.hpp"


namespace EvilPack {


template<class T>
EvilBase<T>::~EvilBase()
{}


template<class T>
RCP<EvilBase<T> >
EvilBase<T>::createEvil(const std::string& concreteEvilName)
{
  if (concreteEvilName == "AEvil") {
    return aEvil<T>();
  }
  else if (concreteEvilName == "BEvil") {
    return bEvil<T>();
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return Teuchos::null; // Never be executed
}


} // namespace EvilPack


#endif // EVIL_BASE_DEF_HPP
