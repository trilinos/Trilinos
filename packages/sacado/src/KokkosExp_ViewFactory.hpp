// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOS_EXPERIMENTAL_VIEW_FACTORY_HPP
#define KOKKOS_EXPERIMENTAL_VIEW_FACTORY_HPP

// This only works with the experimental view enabled
#include "Sacado_ConfigDefs.h"
#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

#include "Sacado_Traits.hpp"
#include "KokkosExp_View_Fad.hpp"

namespace Kokkos {

namespace Impl {

// Class to determine the value_type for a view as a function of one or more
// input views
template <class ... ViewPack>
struct ViewFactoryType {};

template <class View>
struct ViewFactoryType<View> {
  typedef typename View::value_type type;
};

template <class View, class ... ViewPack>
struct ViewFactoryType<View,ViewPack...> {
  typedef typename Sacado::Promote<
    typename View::value_type,
    typename ViewFactoryType<ViewPack...>::type
    >::type type;
};

}

// Function to compute the scalar dimension (e.g., Fad dimesion) from one or
// more views.  It relies on the overload for a single view provided by Sacado
template <class View, class ... ViewPack>
unsigned dimension_scalar(const View& v, const ViewPack&... views) {
  const unsigned dim0 = dimension_scalar(v);
  const unsigned dim1 = dimension_scalar(views...);
  return dim0 >= dim1 ? dim0 : dim1 ;
}

// Traits class used to create a view for a given rank and dimension as a
// function of one or more views.  The value_type for the view is determined
// by value_type, and the view is created through the create_view() function.
// The calling code must determine the rank and dimensions of the view to create
// however internal Sacado dimension will be determined automatically.
template <class ... ViewPack>
struct ViewFactory {

  typedef typename Impl::ViewFactoryType<ViewPack...>::type value_type;

  template <class ResultView, class CtorProp, class ... Dims>
  static ResultView
  create_view(const ViewPack& ... views,
              const CtorProp& prop,
              const Dims ... dims)
    {
      typename ResultView::array_layout layout(dims...);
      layout.dimension[ unsigned(ResultView::rank) ] =
        dimension_scalar(views...);
      return ResultView(prop, layout);
    }
};

}

#endif

#endif /* #ifndef KOKKOS_EXPERIMENTAL_VIEW_FACTORY_HPP */
