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

#ifndef KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP
#define KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP

#include "Kokkos_DynRankView.hpp"

#include "Sacado_ConfigDefs.h"

#if defined(HAVE_SACADO_KOKKOSCONTAINERS) && defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Kokkos_View_Fad.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <>
struct DynRankDimTraits<ViewSpecializeSacadoFad> {

  // Compute the rank of the view from the nonzero dimension arguments.
  // For views of Fad, the rank is one less than the rank determined by the nonzero dimension args
  static size_t computeRank( const size_t N0
                           , const size_t N1
                           , const size_t N2
                           , const size_t N3
                           , const size_t N4
                           , const size_t N5
                           , const size_t N6
                           , const size_t N7 )
  {
    return
      (   (N7 == 0 && N6 == 0 && N5 == 0 && N4 == 0 && N3 == 0 && N2 == 0 && N1 == 0 && N0 == 0) ? 0
      : ( (N7 == 0 && N6 == 0 && N5 == 0 && N4 == 0 && N3 == 0 && N2 == 0 && N1 == 0) ? 0
      : ( (N7 == 0 && N6 == 0 && N5 == 0 && N4 == 0 && N3 == 0 && N2 == 0) ? 1
      : ( (N7 == 0 && N6 == 0 && N5 == 0 && N4 == 0 && N3 == 0) ? 2
      : ( (N7 == 0 && N6 == 0 && N5 == 0 && N4 == 0) ? 3
      : ( (N7 == 0 && N6 == 0 && N5 == 0) ? 4
      : ( (N7 == 0 && N6 == 0) ? 5
      : ( (N7 == 0) ? 6
      : 7 ) ) ) ) ) ) ) );
  }

  // Compute the rank of the view from the nonzero layout arguments.
  template <typename Layout>
  static size_t computeRank( const Layout& layout )
  {
    return computeRank( layout.dimension[0]
                      , layout.dimension[1]
                      , layout.dimension[2]
                      , layout.dimension[3]
                      , layout.dimension[4]
                      , layout.dimension[5]
                      , layout.dimension[6]
                      , layout.dimension[7] );
  }

  // Create the layout for the rank-7 view.
  // For Fad we have to move the fad dimension to the last (rank 8 since the DynRankView is rank-7)
  template <typename Layout>
  static Layout createLayout( const Layout& layout )
  {
    Layout l( layout.dimension[0] != 0 ? layout.dimension[0] : 1
            , layout.dimension[1] != 0 ? layout.dimension[1] : 1
            , layout.dimension[2] != 0 ? layout.dimension[2] : 1
            , layout.dimension[3] != 0 ? layout.dimension[3] : 1
            , layout.dimension[4] != 0 ? layout.dimension[4] : 1
            , layout.dimension[5] != 0 ? layout.dimension[5] : 1
            , layout.dimension[6] != 0 ? layout.dimension[6] : 1
            , layout.dimension[7] != 0 ? layout.dimension[7] : 1 );
    const unsigned fad_dim = computeRank(layout);
    const size_t fad_size = layout.dimension[fad_dim];
    l.dimension[fad_dim] = 1;
    l.dimension[7] = fad_size;

    return l;
  }

  // Create a view from the given dimension arguments.
  // This is only necessary because the shmem constructor doesn't take a layout.
  template <typename ViewType, typename ViewArg>
  static ViewType createView( const ViewArg& arg
                            , const size_t N0
                            , const size_t N1
                            , const size_t N2
                            , const size_t N3
                            , const size_t N4
                            , const size_t N5
                            , const size_t N6
                            , const size_t N7 )
  {
    typename ViewType::array_layout l( N0, N1, N2, N3, N4, N5, N6, N7 );
    typename ViewType::array_layout l_fad = createLayout(l);
    return ViewType( arg
                   , l_fad.dimension[0]
                   , l_fad.dimension[1]
                   , l_fad.dimension[2]
                   , l_fad.dimension[3]
                   , l_fad.dimension[4]
                   , l_fad.dimension[5]
                   , l_fad.dimension[6]
                   , l_fad.dimension[7] );
  }

};

}
}
}

#endif

#include "Kokkos_DynRankView.hpp"

#endif /* #ifndef KOKKOS_DYN_RANK_VIEW_SACADO_FAD_HPP */
