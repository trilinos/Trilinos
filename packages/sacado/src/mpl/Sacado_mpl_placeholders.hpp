// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_PLACEHOLDERS_HPP
#define SACADO_MPL_PLACEHOLDERS_HPP

#include "Sacado_mpl_none.hpp"

namespace Sacado {

  namespace mpl {

    // Placeholder definitions
    template <int N> struct arg {};

    template <> struct arg<1> {
      template <class A1, 
                class A2=mpl::none, 
                class A3=mpl::none, 
                class A4=mpl::none, 
                class A5=mpl::none>
      struct apply {
        typedef A1 type;
      };
    };
    template <> struct arg<2> {
      template <class A1, 
                class A2, 
                class A3=mpl::none, 
                class A4=mpl::none, 
                class A5=mpl::none>
      struct apply {
        typedef A2 type;
      };
    };
    template <> struct arg<3> {
      template <class A1, 
                class A2, 
                class A3, 
                class A4=mpl::none, 
                class A5=mpl::none>
      struct apply {
        typedef A3 type;
      };
    };
    template <> struct arg<4> {
      template <class A1, 
                class A2, 
                class A3, 
                class A4, 
                class A5=mpl::none>
      struct apply {
        typedef A4 type;
      };
    };
    template <> struct arg<5> {
      template <class A1, 
                class A2, 
                class A3, 
                class A4, 
                class A5>
      struct apply {
        typedef A5 type;
      };
    };

    // Placeholder synonyms
    namespace placeholders {
      typedef arg<1> _1;
      typedef arg<2> _2;
      typedef arg<3> _3;
      typedef arg<4> _4;
      typedef arg<5> _5;
      typedef arg<-1> _;
    } // namespace placeholders

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_PLACEHOLDERS_HPP
