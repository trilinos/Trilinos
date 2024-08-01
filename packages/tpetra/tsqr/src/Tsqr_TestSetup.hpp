// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_TestSetup_hpp
#define __TSQR_TestSetup_hpp

#include "Tsqr_MessengerBase.hpp"
#include "Tsqr_Random_GlobalMatrix.hpp"
#include "Tsqr_Matrix.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <vector>

namespace TSQR {
  namespace Test {
    template<class Ordinal, class CommOrdinal>
    Ordinal
    numLocalRows (const Ordinal nrowsGlobal,
                  const CommOrdinal myRank,
                  const CommOrdinal nprocs)
    {
      const Ordinal nrowsLocal = nrowsGlobal / Ordinal(nprocs);
      const Ordinal remainder = nrowsGlobal - nrowsLocal * Ordinal(nprocs);
      if (myRank != nprocs - 1)
        return nrowsLocal;
      else
        return nrowsLocal + remainder;
    }

    /// \param generator [in/out] Proc 0 is the only MPI process that
    /// generates pseudorandom numbers.  This allows us to use a
    /// sequential PRNG.  Otherwise, we would have to use a parallel
    /// PRNG, in order to prevent correlations between numbers
    /// generated on different MPI processes.  On processes other than
    /// Proc 0, generator is not touched.
    template<class MatrixViewType, class Generator>
    void
    distributedTestProblem (Generator& generator,
                            MatrixViewType& A_local,
                            MessengerBase<typename MatrixViewType::ordinal_type>* const ordinalComm,
                            MessengerBase<typename MatrixViewType::non_const_value_type>* const scalarComm)
    {
      using ordinal_type = typename MatrixViewType::ordinal_type;
      using scalar_type =
        typename MatrixViewType::non_const_value_type;
      using magnitude_type =
        typename Teuchos::ScalarTraits<scalar_type>::magnitudeType;

      const int myRank = scalarComm->rank();
      const ordinal_type ncols = A_local.extent(1);
      if (myRank == 0) {
        // Generate some singular values for the test problem.
        std::vector<magnitude_type> singular_values (ncols);
        singular_values[0] = 1.0;
        for (ordinal_type k = 1; k < ncols; ++k) {
          singular_values[k] = singular_values[k-1] / magnitude_type(2.0);
        }

        // Generate the test problem.  All MPI processes
        // participate, but only Proc 0 generates the (pseudo)random
        // numbers.
        TSQR::Random::randomGlobalMatrix (&generator, A_local,
                                          singular_values.data (),
                                          ordinalComm, scalarComm);
      }
      else {
        // This helps C++ deduce the type; the values aren't read on
        // this proc.
        magnitude_type singular_values[1];

        // All MPI processes participate in the distribution of the
        // test matrix.
        TSQR::Random::randomGlobalMatrix (&generator, A_local,
                                          singular_values,
                                          ordinalComm, scalarComm);
      }
    }
  } // namespace Test
} // namespace TSQR

#endif // __TSQR_TestSetup_hpp
