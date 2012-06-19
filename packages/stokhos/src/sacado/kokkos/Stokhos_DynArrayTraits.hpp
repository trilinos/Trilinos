// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_DYN_ARRAY_TRAITS_HPP
#define STOKHOS_DYN_ARRAY_TRAITS_HPP

namespace Stokhos {

  //! Base template specification for %IsScalarType
  /*!
   * The %IsScalarType classes provide a mechanism for computing the 
   * determining whether a type is a scalar type (float, double, etc...)
   */
  template <typename T> struct IsScalarType2 {
    static const bool value = false;
  };

  //! Specialization of above classes to built-in types
#define STOKHOS_BUILTIN_SPECIALIZATION(t)                  \
  template <> struct IsScalarType2< t > {	          \
    static const bool value = true;	       		  \
  };

  STOKHOS_BUILTIN_SPECIALIZATION(float)
  STOKHOS_BUILTIN_SPECIALIZATION(double)
  STOKHOS_BUILTIN_SPECIALIZATION(int)
  STOKHOS_BUILTIN_SPECIALIZATION(long)

#undef STOKHOS_BUILTIN_SPECIALIZATION
  

  /*!
   * \brief Dynamic array allocation class that works for any type
   */
  template <typename T, typename node, bool isScalar = IsScalarType2<T>::value>
  struct DynArrayTraits {

    //! Copy array from \c src to \c dest of length \c sz
    static inline void copy(const T* src, T*  dest, std::size_t sz);

    //! Zero out array \c dest of length \c sz
    static inline void zero(T* dest, std::size_t sz);

    //! Fill array \c dest of length \c sz with value \c v
    static inline void fill(T* desk, std::size_t sz, const T& v);

    //! Get memory for new array of length \c sz and fill with zeros
    static inline T* get_and_fill(std::size_t sz);

    /*! 
     * \brief Get memory for new array of length \c sz and fill with 
     * entries from \c src
     */
    static inline T* get_and_fill(const T* src, std::size_t sz);

    //! Destroy array elements and release memory
    static inline void destroy_and_release(T* m, std::size_t sz);
  };

} // namespace Stokhos

// Host specialization
#include "KokkosArray_Host.hpp"
#include "KokkosArray_Host_macros.hpp"
#include "Stokhos_DynArrayTraits_impl.hpp"
#include "KokkosArray_Clear_macros.hpp"

// Cuda specialization
#include "KokkosArray_Cuda.hpp"
#include "KokkosArray_Cuda_macros.hpp"
#include "Stokhos_DynArrayTraits_impl.hpp"
#include "KokkosArray_Clear_macros.hpp"

#endif // STOKHOS_DYN_ARRAY_TRAITS_HPP
