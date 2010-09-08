/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

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
