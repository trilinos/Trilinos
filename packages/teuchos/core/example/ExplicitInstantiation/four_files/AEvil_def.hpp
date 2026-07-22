// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef A_EVIL_DEF_HPP
#define A_EVIL_DEF_HPP


#include "AEvil_decl.hpp"
// We have to include this file to be 100% safe since we included
// EvilBaseDecl.hpp in AEvilDecl.hpp
#include "EvilBase.hpp"
// We need to have BEvil's implementation to call it!
#include "BEvil.hpp"


namespace EvilPack {


template<class T>
void AEvil<T>::callBEvil(const BEvil<T> &bEvil, const T& obj) const
{
  using Teuchos::typeName;
  std::cout << typeName(*this) << " call BEvil: ";
  bEvil.soundOff(obj);
}


template<class T>
void AEvil<T>::soundOff(const T& obj) const
{
  using Teuchos::typeName;
  std::cout << typeName(*this) << " obj = " << obj << "\n";
}


} // namespace EvilPack


#endif // A_EVIL_DEF_HPP
