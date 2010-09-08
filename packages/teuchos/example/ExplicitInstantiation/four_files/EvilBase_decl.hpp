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

#ifndef EVIL_BASE_DECL_HPP
#define EVIL_BASE_DECL_HPP


#include "Teuchos_RCP.hpp"


namespace EvilPack {


/** \brief . */
using Teuchos::RCP;


/** \brief Evil base class that people often write with a factory function
 * to all of the subclasses in the interface.
 */
template<class T>
class EvilBase {
public:

  /** \brief. Required virtual destructor. */
  virtual ~EvilBase();

  /** \brief The virtual function. */
  virtual void soundOff(const T& obj) const = 0;
  
  /** \brief The factory in the interface. */
  static RCP<EvilBase<T> >
  createEvil(const std::string& concreteEvilName);

};


} // namespace EvilPack


#endif // EVIL_BASE_DECL_HPP
