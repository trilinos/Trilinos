// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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
    SACADO_INLINE_FUNCTION
    static void copy(const T* src, T*  dest, int sz) {
      for (int i=0; i<sz; ++i)
        *(dest++) = *(src++);
    }

    //! Zero out array \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void zero(T* dest, int sz) {
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
    SACADO_INLINE_FUNCTION
    static void copy(const T* src, T* dest, int sz) {
      if (sz > 0)
#if defined( __CUDACC__) || defined(__HIPCC__ )
        for (int i=0; i<sz; ++i)
          dest[i] = src[i];
#else
        std::memcpy(dest,src,sz*sizeof(T));
#endif
    }

    //! Zero out array \c dest of length \c sz
    SACADO_INLINE_FUNCTION
    static void zero(T* dest, int sz) {
      if (sz > 0)
#if defined(__CUDACC__ ) || defined(__HIPCC__ )
        for (int i=0; i<sz; ++i)
          dest[i] = T(0.);
#else
        std::memset(dest,0,sz*sizeof(T));
#endif
    }

  };

} // namespace Sacado

#endif // SACAD0_STATICARRAYTRAITS_HPP
