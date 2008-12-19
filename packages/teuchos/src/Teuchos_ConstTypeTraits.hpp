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

#ifndef TEUCHOS_CONST_TYPE_TRAITS_HPP
#define TEUCHOS_CONST_TYPE_TRAITS_HPP


#include "Teuchos_ConfigDefs.hpp"


namespace Teuchos {


/** \brief Traits class that strips 'const' off of a type.
 *
 * See "Modern C++ Design", 2001, by Andrei Alexandrescu.
 *
 * \ingroup teuchos_language_support_grp
 */
template<class T>
class ConstTypeTraits {
private:
  /** \brief . */
  template<class U> struct UnConst
  { typedef U Result; };
  /** \brief . */
  template<class U> struct UnConst<const U>
  { typedef U Result; };
public:
  /** \brief . */
  typedef typename UnConst<T>::Result NonConstType;
};


} // end namespace Teuchos


#endif	// TEUCHOS_CONST_TYPE_TRAITS_HPP
