// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef B_EVIL_DEF_HPP
#define B_EVIL_DEF_HPP


#include "BEvil_decl.hpp"
// We have to include this file to be 100% safe since we included
// EvilBaseDecl.hpp in BEvilDecl.hpp
#include "EvilBase.hpp"
// We need to have AEvil's implementation to call it!
#include "AEvil.hpp"


namespace EvilPack {


template<class T>
void BEvil<T>::callAEvil(const AEvil<T> &aEvil, const T& obj) const
{
  using Teuchos::typeName;
  std::cout << typeName(*this) << " call AEvil: ";
  aEvil.soundOff(obj);
}


template<class T>
void BEvil<T>::soundOff(const T& obj) const
{
  using Teuchos::typeName;
  std::cout << typeName(*this) << " obj = " << obj << "\n";
}


} // namespace EvilPack


#endif // B_EVIL_DEF_HPP
