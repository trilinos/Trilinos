// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_GETRAWPTR_HPP
#define TEUCHOS_GETRAWPTR_HPP

#include "Teuchos_ArrayRCP.hpp"

namespace Teuchos {

template<class Container>
class RawPointerConversionTraits {
  // Empty : will fail if specialization doesn't exist
};


template<class Container>
typename RawPointerConversionTraits<Container>::Ptr_t
getRawPtr( const Container& c )
{
  return RawPointerConversionTraits<Container>::getRawPtr(c);
}


// partial specialization for C pointer
template<class RawType>
class RawPointerConversionTraits<RawType*>
{
public:
  typedef RawType* Ptr_t;
  static Ptr_t getRawPtr( RawType* p ) { return p; }
};

// partial specialization for ArrayRCP
template<class T>
class RawPointerConversionTraits<ArrayRCP<T> >
{
public:
  typedef typename ArrayRCP<T>::pointer Ptr_t;
  static Ptr_t getRawPtr( const ArrayRCP<T>& arcp ) { return arcp.getRawPtr(); }
};

// ToDo: Add specializations as needed!

} // namespace Teuchos

#endif // TEUCHOS_GETRAWPTR_HPP

