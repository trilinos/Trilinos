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

#ifndef _TEUCHOS_TYPE_TRAITS_HPP_
#define _TEUCHOS_TYPE_TRAITS_HPP_

/*! \file Teuchos_TypeTraits.hpp
 \brief Defines basic traits allowing evaluation of type equality.
*/

namespace  Teuchos {

namespace TypeTraits {

/** \brief Default is_equal traits class has \c value equal to \c false, indicating that T1 and T2 are not equal, 
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T1, typename T2>
struct is_same {
  enum {value = false};
};

/** \brief Partial specialization of is_equal class for equal types, where \c value equal to \c true.
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T>
struct is_same<T,T> {
  enum {value = true};
};
 
}

} // namespace Teuchos


#endif // _TEUCHOS_TYPE_TRAITS_HPP_

