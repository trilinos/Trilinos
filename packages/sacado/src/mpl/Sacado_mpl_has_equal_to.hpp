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

#ifndef SACADO_MPL_HAS_EQUAL_TO_HPP
#define SACADO_MPL_HAS_EQUAL_TO_HPP

#include <type_traits>

#include "Sacado_mpl_void.hpp"

namespace Sacado {

  namespace mpl {

    template <typename T1, typename T2 = T1, typename = void_t<> >
    struct has_equal_to : std::false_type {};

    template <typename T1, typename T2>
    struct has_equal_to<T1, T2, void_t<decltype(std::declval<T1>() ==
                                                std::declval<T2>())> >
    : std::true_type {};

  }

}

#endif // SACADO_MPL_HAS_EQUAL_TO_HPP
