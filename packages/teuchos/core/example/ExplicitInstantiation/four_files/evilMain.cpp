// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


// The client just includes the basic *.hpp forms without having the worry
// about implicit or explicit instantiation!
#include "EvilBase.hpp"
#include "AEvil.hpp"
#include "BEvil.hpp"


template<class T>
void testEvil(const T& obj)
{

  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using EvilPack::EvilBase;
  using EvilPack::AEvil;
  using EvilPack::BEvil;

  RCP<AEvil<T> > aEvil =
    rcp_dynamic_cast<AEvil<T> >(EvilBase<T>::createEvil("AEvil"));
  RCP<BEvil<T> > bEvil =
    rcp_dynamic_cast<BEvil<T> >(EvilBase<T>::createEvil("BEvil"));

  aEvil->soundOff(obj);
  bEvil->soundOff(obj);
  aEvil->callBEvil(*bEvil, obj);
  bEvil->callAEvil(*aEvil, obj);

}


int main()
{
  testEvil<double>(1.0);
  testEvil<int>(2);
  return 0;
}
