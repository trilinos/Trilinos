// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_STATICARRAYTRAITS_HPP
#define SACADO_STATICARRAYTRAITS_HPP

#include <cstring>

#include "Sacado_Traits.hpp"

namespace Sacado {

  /*!
   * \brief Static array allocation class that works for any type
   */
  template <typename T, bool isScalar = IsScalarType<T>::value>
  struct ss_array {
    
    //! Copy array from \c src to \c dest of length \c sz
    static inline void copy(const T* src, T*  dest, int sz) {
      for (int i=0; i<sz; ++i)
	*(dest++) = *(src++);
    }

    //! Zero out array \c dest of length \c sz
    static inline void zero(T* dest, int sz) {
      for (int i=0; i<sz; ++i)
	*(dest++) = T(0.);
    }
  };

  /*!
   * \brief Static array allocation class that is specialized for scalar
   * i.e., fundamental or built-in types (float, double, etc...).
   */
  template <typename T>
  struct ss_array<T,true> {
    
    //! Copy array from \c src to \c dest of length \c sz
    static inline void copy(const T* src, T* dest, int sz) {
      std::memcpy(dest,src,sz*sizeof(T));
    }
    
    //! Zero out array \c dest of length \c sz
    static inline void zero(T* dest, int sz) {
      std::memset(dest,0,sz*sizeof(T));
    }

  };

} // namespace Sacado

#endif // SACAD0_STATICARRAYTRAITS_HPP
