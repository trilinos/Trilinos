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

#ifndef SCADO_MPL_VECTOR_AT_SPEC_HPP
#define SCADO_MPL_VECTOR_AT_SPEC_HPP

namespace Sacado {

  namespace mpl {

    template <class Vector, int Pos> struct vector_at {};

    template <typename T, typename...Args>
    struct vector_at<mpl::vector<T,Args...>, 0> {
      typedef T type;
    };

    template <typename T, typename...Args, int Pos>
    struct vector_at<mpl::vector<T,Args...>, Pos> {
      typedef typename vector_at<mpl::vector<Args...>, Pos-1>::type type;
    };

  }

}

#endif // SCADO_MPL_VECTOR_AT_SPEC_HPP
