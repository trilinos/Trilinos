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

#ifndef SACADO_MPL_APPLY_HPP
#define SACADO_MPL_APPLY_HPP

#include "Sacado_mpl_apply_wrap.hpp"
#include "Sacado_mpl_lambda.hpp"
#include "Sacado_mpl_none.hpp"

namespace Sacado {

  namespace mpl {

    template <class F> 
    struct apply0 : apply_wrap0<typename lambda<F>::type> {};

    template <class F, class A1> 
    struct apply1 : apply_wrap1<typename lambda<F>::type,A1> {};

    template <class F, class A1, class A2> 
    struct apply2 : apply_wrap2<typename lambda<F>::type,A1,A2> {};

    template <class F, class A1, class A2, class A3> 
    struct apply3 : apply_wrap3<typename lambda<F>::type,A1,A2,A3> {};

    template <class F, class A1, class A2, class A3, class A4> 
    struct apply4 : apply_wrap4<typename lambda<F>::type,A1,A2,A3,A4> {};

    template <class F, class A1, class A2, class A3, class A4, class A5> 
    struct apply5 : apply_wrap5<typename lambda<F>::type,A1,A2,A3,A4,A5> {};

    template <class F, 
              class A1=mpl::none, 
              class A2=mpl::none, 
              class A3=mpl::none, 
              class A4=mpl::none, 
              class A5=mpl::none>
    struct apply : apply_wrap<typename lambda<F>::type,A1,A2,A3,A4,A5> {};

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_APPLY_HPP
