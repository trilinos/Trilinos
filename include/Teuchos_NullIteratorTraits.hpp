// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
