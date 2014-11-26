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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_MPL_TYPE_WRAP_HPP
#define SACADO_MPL_TYPE_WRAP_HPP

#include "Sacado_mpl_has_type.hpp"

namespace Sacado {

  namespace mpl {

    // Wrap a type T so that it becomes a metafunction if it isn't already

    // Don't use mpl_if because it uses type_wrap in its implementation
    template <class T> struct add_type { typedef T type; };
    template <bool cond, class T> struct type_wrap_impl {};
    template <class T> struct type_wrap_impl<true,T> : T {};
    template <class T> struct type_wrap_impl<false,T> : add_type<T> {};
    template <class T> struct type_wrap :
      type_wrap_impl< mpl::has_type<T>::value, T > {};

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_TYPE_WRAP_HPP
