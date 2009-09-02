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

