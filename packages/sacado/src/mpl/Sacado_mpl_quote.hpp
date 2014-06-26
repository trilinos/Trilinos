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

#ifndef SACADO_MPL_QUOTE_HPP
#define SACADO_MPL_QUOTE_HPP

#include "Sacado_mpl_type_wrap.hpp"

namespace Sacado {

  namespace mpl {
   
    // Transform a class/template to a metafunction class 
    // (nested template apply)
    template <class F>
    struct quote0 {
      struct apply : mpl::type_wrap< F > {};
    };

    template <template<class T1> class F>
    struct quote1 {
      template <class T1> 
      struct apply : mpl::type_wrap< F<T1> > {};
    };

    template < template<class T1, 
                        class T2> class F>
    struct quote2 {
      template <class T1, 
                class T2> 
      struct apply : mpl::type_wrap< F<T1,T2> >{};
    };

    template < template<class T1, 
                        class T2, 
                        class T3> class F>
    struct quote3 {
      template <class T1, 
                class T2, 
                class T3> 
      struct apply : mpl::type_wrap< F<T1,T2,T3> >{};
    };

    template < template<class T1, 
                        class T2, 
                        class T3, 
                        class T4> class F>
    struct quote4 {
      template <class T1, 
                class T2, 
                class T3,
                class T4> 
      struct apply : mpl::type_wrap< F<T1,T2,T3,T4> >{};
    };

    template < template<class T1, 
                        class T2, 
                        class T3, 
                        class T4, 
                        class T5> class F>
    struct quote5 {
      template <class T1, 
                class T2, 
                class T3,
                class T4,
                class T5> 
      struct apply : mpl::type_wrap< F<T1,T2,T3,T4,T5> >{};
    };

    template <class F>
    struct quote : quote0<F> {};

    template < template<class T1> class F,
               class T1>
    struct quote< F<T1> > : quote1<F> {};

    template < template<class T1, 
                        class T2> class F,
               class T1,
               class T2>
    struct quote< F<T1,T2> > : quote2<F> {};

    template < template<class T1, 
                        class T2, 
                        class T3> class F,
               class T1,
               class T2,
               class T3>
    struct quote< F<T1,T2,T3> > : quote3<F> {};

    template < template<class T1, 
                        class T2, 
                        class T3, 
                        class T4> class F,
               class T1,
               class T2,
               class T3,
               class T4>
    struct quote< F<T1,T2,T3,T4> > : quote4<F> {};

    template < template<class T1, 
                        class T2, 
                        class T3, 
                        class T4, 
                        class T5> class F,
               class T1,
               class T2,
               class T3,
               class T4,
               class T5>
    struct quote< F<T1,T2,T3,T4,T5> > : quote5<F> {};

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_QUOTE_HPP
