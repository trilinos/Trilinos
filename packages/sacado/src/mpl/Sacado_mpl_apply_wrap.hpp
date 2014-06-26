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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_MPL_APPLY_WRAP_HPP
#define SACADO_MPL_APPLY_WRAP_HPP

#include "Sacado_mpl_none.hpp"

namespace Sacado {

  namespace mpl {

    // wrapper to call a metafunction class
    template <class F> struct 
    apply_wrap0 {
      typedef typename F::apply::type type;
    };

    template <class F, 
              class A1> 
    struct apply_wrap1 {
      typedef typename F::template apply<A1>::type type;
    };

    template <class F, 
              class A1, 
              class A2> 
    struct apply_wrap2 {
      typedef typename F::template apply<A1,A2>::type type;
    };

    template <class F, 
              class A1, 
              class A2, 
              class A3> 
    struct apply_wrap3 {
      typedef typename F::template apply<A1,A2,A3>::type type;
    };

    template <class F, 
              class A1, 
              class A2, 
              class A3,
              class A4> 
    struct apply_wrap4 {
      typedef typename F::template apply<A1,A2,A3,A4>::type type;
    };

    template <class F, 
              class A1, 
              class A2, 
              class A3,
              class A4,
              class A5> 
    struct apply_wrap5 {
      typedef typename F::template apply<A1,A2,A3,A4,A5>::type type;
    };

    template <class F, 
              class A1, 
              class A2, 
              class A3,
              class A4,
              class A5>
    struct apply_wrap : 
      apply_wrap5<F,A1,A2,A3,A4,A5> {};

    template <class F, 
              class A1, 
              class A2, 
              class A3,
              class A4>
    struct apply_wrap<F,A1,A2,A3,A4,mpl::none> : 
      apply_wrap4<F,A1,A2,A3,A4> {};

    template <class F, 
              class A1, 
              class A2, 
              class A3>
    struct apply_wrap<F,A1,A2,A3,mpl::none,mpl::none> : 
      apply_wrap3<F,A1,A2,A3> {};

    template <class F, 
              class A1, 
              class A2>
    struct apply_wrap<F,A1,A2,mpl::none,mpl::none,mpl::none> : 
      apply_wrap2<F,A1,A2> {};

    template <class F, 
              class A1>
    struct apply_wrap<F,A1,mpl::none,mpl::none,mpl::none,mpl::none> : 
      apply_wrap1<F,A1> {};

    template <class F>
    struct apply_wrap<F,mpl::none,mpl::none,mpl::none,mpl::none,mpl::none> : 
      apply_wrap0<F> {};
    
  } // namespace mpl

} // namespace Sacado

#endif // SACADO_APPLY_TYPE_WRAP_HPP
