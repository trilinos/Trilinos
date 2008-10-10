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

#ifndef TEUCHOS_NULL_ITERATOR_TRAITS_HPP
#define TEUCHOS_NULL_ITERATOR_TRAITS_HPP


#include "Teuchos_ConfigDefs.hpp"


namespace Teuchos {



/** \brief Base traits class for getting a properly initialized null pointer.
 *
 * This default traits class simply passes in a raw '0' to the constructor.
 * This works for a amazingly large number of classes that can be used a an
 * iterator.
 *
 * \ingroup teuchos_mem_mng_grp
 */
template<typename Iter>
class NullIteratorTraits {
public:
  static Iter getNull()
  {
#ifdef TEUCHOS_NO_ZERO_ITERATOR_CONVERSION
    return Iter();
#else
    return Iter(0);
#endif
  }
};


/** \brief Partial specialization for std::reverse_iterator. 
 *
 * \ingroup teuchos_mem_mng_grp
 */
template<typename Iter>
class NullIteratorTraits<std::reverse_iterator<Iter> > {
public:
  static std::reverse_iterator<Iter> getNull()
    {
      return std::reverse_iterator<Iter>(
        NullIteratorTraits<Iter>::getNull()
        );
    }
};


} // namespace Teuchos


#endif // TEUCHOS_NULL_ITERATOR_TRAITS_HPP
