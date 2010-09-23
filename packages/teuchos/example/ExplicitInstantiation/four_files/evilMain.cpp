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
