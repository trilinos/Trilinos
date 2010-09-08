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

#ifndef B_EVIL_DECL_HPP
#define B_EVIL_DECL_HPP

// Only include the declaration, not any implementations in case of cicular
// dependencies!
#include "EvilBase_decl.hpp"


namespace EvilPack {


// Need a forward for B to declare function callAEvil(...)
template<class T> class AEvil;


/** \brief A subclass of EvilBase that calls AEvil.
 */
template<class T>
class BEvil : public EvilBase<T> {
public:
  /** \brief . */
  void callAEvil(const AEvil<T> &aEvil, const T& obj) const;
  /** \brief . */
  void soundOff(const T& obj) const;
};


/** \brief Nonmember constructor.
 *
 * \relates BEvil
 **/
template<class T>
inline
RCP<BEvil<T> > bEvil()
{
  return Teuchos::rcp(new BEvil<T>);
}


} // namespace EvilPack


#endif // B_EVIL_DECL_HPP
