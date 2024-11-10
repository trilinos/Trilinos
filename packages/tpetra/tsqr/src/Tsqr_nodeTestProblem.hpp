// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Test_nodeTestProblem_hpp
#define __TSQR_Test_nodeTestProblem_hpp

#include "Tsqr_Random_MatrixGenerator.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <algorithm>
#include <vector>

namespace TSQR {
  namespace Test {

    /// Fill in the nrows by ncols matrix A (with leading dimension
    /// lda) with a test problem for single-node TSQR.
    template< class Ordinal, class Scalar, class Generator >
    void
    nodeTestProblem (Generator& generator,
                     const Ordinal nrows,
                     const Ordinal ncols,
                     Scalar A[],
                     const Ordinal lda,
                     const bool numerically_interesting)
    {
      typedef TSQR::Random::MatrixGenerator<Ordinal, Scalar, Generator> matgen_type;
      typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

      matgen_type matGen (generator);

      if (numerically_interesting) {
        std::vector<magnitude_type> singular_values (std::min (nrows, ncols));
        singular_values[0] = magnitude_type (1);
        for (Ordinal k = 1; k < std::min(nrows, ncols); ++k) {
          singular_values[k] = singular_values[k-1] / magnitude_type(2);
        }
        matGen.fill_random_svd (nrows, ncols, A, lda, &singular_values[0]);
      }
      else {
        matGen.fill_random (nrows, ncols, A, lda);
      }
    }

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_nodeTestProblem_hpp
