// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_PUSH_BACK_HPP
#define SACADO_MPL_PUSH_BACK_HPP

namespace Sacado {

  namespace mpl {

    template <class T> struct push_back_impl {};

    template <class Seq, class T>
    struct push_back : 
      push_back_impl<typename Seq::tag>:: template apply<Seq,T> {};

  }

}

#endif // SACADO_MPL_PUSH_BACK_HPP
