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

#ifndef SACADO_MPL_IS_CONVERTIBLE_HPP
#define SACADO_MPL_IS_CONVERTIBLE_HPP

namespace Sacado {

  namespace mpl {

    //
    // A simplified implementation of boost type-trait
    // is_convertible<From,To>.  We use this in a much more limited context
    // within Sacado, and so the easy implementation should always work.
    // We assume From and To are "scalar" types, e.g., are not pointer or
    // reference types.
    //

    struct convertible_impl {
      typedef char yes;       // sizeof(yes) == 1
      typedef char (&no)[2];  // sizeof(no)  == 2

      // A function that takes anything convertible to a To
      template <typename To> static yes tester(To);

      // Overload resolution prefers anything over ...
      template <typename To> static no tester(...);

      // Check if From is convertible to To
      template <typename From, typename To>
      struct checker {
        static From& f;
        static const bool value = sizeof(tester<To>(f)) == sizeof(yes);
      };
    };

    template <typename From, typename To>
    struct is_convertible {
      static const bool value = convertible_impl::checker<From,To>::value;
    };

  }

}

#endif // SACADO_MPL_IS_CONVERTIBLE_HPP
