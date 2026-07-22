// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_RANGE_C_HPP
#define SACADO_MPL_RANGE_C_HPP

#include "Sacado_mpl_none.hpp"
#include "Sacado_mpl_size.hpp"
#include "Sacado_mpl_begin.hpp"
#include "Sacado_mpl_end.hpp"
#include "Sacado_mpl_next.hpp"
#include "Sacado_mpl_at.hpp"
#include "Sacado_mpl_deref.hpp"
#include "Sacado_mpl_integral_c.hpp"

namespace Sacado {

  namespace mpl {

    // range_c tag for mpl operations
    struct range_c_tag {};

    // range_c
    template <class T, T N, T M, T Delta = 1>
    struct range_c {
      typedef range_c_tag tag;
      typedef range_c type;
      static const int sz = (M-N+Delta-1)/Delta;
      typedef T integral_type;
      static const int start_value = N;
      static const int end_value = M;
      static const int step_value = Delta;
    };

    // iterator
    template <class Range, int Pos>
    struct range_c_iterator {
      static const int value = Pos;
    };

    // size
    template <>
    struct size_impl<range_c_tag> {
      template <class Range>
      struct apply {
        static const int value = Range::sz;
      };
    };

    // begin
    template <>
    struct begin_impl<range_c_tag> {
      template <class Range>
      struct apply {
        typedef range_c_iterator<Range,0> type;
      };
    };

    // end
    template <>
    struct end_impl<range_c_tag> {
      template <class Range>
      struct apply {
        typedef range_c_iterator<Range,Range::sz> type;
      };
    };

    // next
    template <class Range, int Pos>
    struct next< range_c_iterator<Range,Pos> > {
      typedef range_c_iterator<Range,Pos+1> type;
    };



    // at
    template <int Pos>
    struct at_impl<range_c_tag, Pos> {
      template <class Range>
      struct apply {
        typedef integral_c<typename Range::integral_type,
                           Range::start_value + Range::step_value*Pos> type;
      };
    };

    // deref
    template <class Range, int Pos>
    struct deref< range_c_iterator<Range,Pos> > : mpl::at<Range,Pos> {};

  }
}

#endif
