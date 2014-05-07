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

#ifndef SACADO_MPL_IS_PLACEHOLDER_HPP
#define SACADO_MPL_IS_PLACEHOLDER_HPP

#include "Sacado_mpl_is_placeholder.hpp"
#include "Sacado_mpl_none.hpp"

namespace Sacado {

  namespace mpl {

    template <class F> 
    struct is_placeholder { 
      static const bool value = false; 
    };
    template <int N> 
    struct is_placeholder< arg<N> > {
      static const bool value = true; 
    };
    template <template <class T1> class F,
              class T1> 
    struct is_placeholder< F<T1> > {
      static const bool value = is_placeholder<T1>::value;
    };
    template <template <class T1, class T2> class F,
              class T1,
              class T2> 
    struct is_placeholder< F<T1,T2> > {
      static const bool value = 
        is_placeholder<T1>::value || 
        is_placeholder<T2>::value;
    };
    template <template <class T1, class T2, class T3> class F,
              class T1,
              class T2,
              class T3> 
    struct is_placeholder< F<T1,T2,T3> > {
      static const bool value = 
        is_placeholder<T1>::value || 
        is_placeholder<T2>::value || 
        is_placeholder<T3>::value;
    };
    template <template <class T1, class T2, class T3, class T4> class F,
              class T1,
              class T2,
              class T3,
              class T4> 
    struct is_placeholder< F<T1,T2,T3,T4> > {
      static const bool value = 
        is_placeholder<T1>::value || 
        is_placeholder<T2>::value || 
        is_placeholder<T3>::value || 
        is_placeholder<T4>::value;
    };
    template <template <class T1, class T2, class T3, class T4, class T5> class F,
              class T1,
              class T2,
              class T3,
              class T4,
              class T5> 
    struct is_placeholder< F<T1,T2,T3,T4,T5> > {
      static const bool value = 
        is_placeholder<T1>::value || 
        is_placeholder<T2>::value || 
        is_placeholder<T3>::value || 
        is_placeholder<T4>::value || 
        is_placeholder<T5>::value;
    };

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_IS_PLACEHOLDER_HPP
