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

#ifndef STOKHOS_STATIC_ARRAY_TRAITS_HPP
#define STOKHOS_STATIC_ARRAY_TRAITS_HPP

#include "Sacado_Traits.hpp"

namespace Stokhos {

  /*!
   * \brief Static array allocation class
   */
  template <typename T, typename node, 
	    bool isScalar = Sacado::IsScalarType<T>::value>
  struct StaticArrayTraits {

    typedef T value_type;
    typedef node node_type;
    
    //! Copy array from \c src to \c dest of length \c sz
    static inline void copy(const T* src, T*  dest, std::size_t sz);

    //! Zero out array \c dest of length \c sz
    static inline void zero(T* dest, std::size_t sz);

    //! Fill array \c dest of length \c sz with value \c v
    static inline void fill(T* desk, std::size_t sz, const T& v);
  };

} // namespace Stokhos

// Host specialization
#include "KokkosArray_Host.hpp"
#include "KokkosArray_Host_macros.hpp"
#include "Stokhos_StaticArrayTraits_impl.hpp"
#include "KokkosArray_Clear_macros.hpp"

// Cuda specialization
#include "KokkosArray_Cuda.hpp"
#include "KokkosArray_Cuda_macros.hpp"
#include "Stokhos_StaticArrayTraits_impl.hpp"
#include "KokkosArray_Clear_macros.hpp"

#endif // STOKHOS_STATIC_ARRAY_TRAITS_HPP
