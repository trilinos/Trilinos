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

#ifndef SACADO_MPL_VECTOR_HPP
#define SACADO_MPL_VECTOR_HPP

#include "Sacado_ConfigDefs.h"
#ifdef HAVE_SACADO_CXX11

#include "Sacado_mpl_none.hpp"
#include "Sacado_mpl_size.hpp"
#include "Sacado_mpl_begin.hpp"
#include "Sacado_mpl_end.hpp"
#include "Sacado_mpl_next.hpp"
#include "Sacado_mpl_push_back.hpp"
#include "Sacado_mpl_at.hpp"
#include "Sacado_mpl_deref.hpp"

// Specializations needed for various mpl operations
#include "Sacado_mpl_vector_size_spec.hpp"

namespace Sacado {

  namespace mpl {

    // vector tag for mpl operations
    struct vector_tag {};

    // vector
    template <typename... Args>
    struct vector :
      vector_size<Args...> {
      typedef vector_tag tag;
      typedef vector type;
    };

    // iterator
    template <class Vector, int Pos>
    struct vector_iterator {
      static const int value = Pos;
    };

    // size
    template <>
    struct size_impl<vector_tag> {
      template <class Vector>
      struct apply {
        static const int value = Vector::sz;
      };
    };

    // begin
    template <>
    struct begin_impl<vector_tag> {
      template <class Vector>
      struct apply {
        typedef vector_iterator<Vector,0> type;
      };
    };

    // end
    template <>
    struct end_impl<vector_tag> {
      template <class Vector>
      struct apply {
        typedef vector_iterator<Vector,Vector::sz> type;
      };
    };

    // next
    template <class Vector, int Pos>
    struct next< vector_iterator<Vector,Pos> > {
      typedef vector_iterator<Vector,Pos+1> type;
    };

    // deref
    template <class Vector, int Pos>
    struct deref< vector_iterator<Vector,Pos> > : mpl::at<Vector,Pos> {};

  }

}

#include "Sacado_mpl_vector_at_spec.hpp"
#include "Sacado_mpl_vector_push_back_spec.hpp"

namespace Sacado {
  namespace mpl {

    // at
    template <int Pos>
    struct at_impl<vector_tag, Pos> {
      template <class Vector>
      struct apply : vector_at<Vector,Pos> {};
    };

    // push_back
    template <>
    struct push_back_impl<vector_tag> {
      template <class Vector, class T>
      struct apply : vector_push_back<Vector, T> {};
    };

  }
}

#endif

#endif
