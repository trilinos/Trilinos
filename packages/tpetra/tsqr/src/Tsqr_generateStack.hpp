// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Test_generateStack_hpp
#define __TSQR_Test_generateStack_hpp

#include "Tsqr_Matrix.hpp"
#include "Tsqr_Util.hpp"
#include "Tsqr_Random_MatrixGenerator.hpp"
#include "Tsqr_RMessenger.hpp"

#include <algorithm>
#include <functional>
#include <sstream>
#include <stdexcept>

namespace TSQR {
  namespace Test {
    /// \brief Generate a random "R stack" test problem on one MPI process.
    ///
    /// Generate a (pseudo)random test problem consisting of numProcs
    /// different numRows by numCols upper triangular matrices,
    /// stacked vertically.
    ///
    /// \param generator [in/out] (Pseudo)random number generator,
    ///   that generates according to a normal(0,1) distribution.
    ///
    /// \param A_global [out] Matrix to fill.  Should be empty on
    ///   input.  Output will be on Proc 0 only.  We will set
    ///   dimensions, allocate, and fill.
    ///
    /// \param singularValues [in] ncols singular values to use
    ///
    /// \param numProcs Number of (MPI) processes
    ///
    /// \param numCols Number of columns in the output matrix A_global
    ///
    template<class Ordinal, class Scalar, class Generator>
    static void
    generateStack (Generator& generator,
                   Matrix<Ordinal, Scalar>& A_global,
                   const typename Teuchos::ScalarTraits<Scalar>::magnitudeType singularValues[],
                   const int numProcs,
                   const Ordinal numCols)
    {
      using mat_view_type = MatView<Ordinal, Scalar>;

      TSQR::Random::MatrixGenerator<Ordinal, Scalar, Generator> matGen (generator);
      const Ordinal numRows = numProcs * numCols;
      A_global.reshape (numRows, numCols);
      deep_copy (A_global, Scalar {});

      for (int p = 0; p < numProcs; ++p) {
        auto* const curptr = A_global.data() + p*numCols;
        mat_view_type R_cur (numCols, numCols, curptr, numRows);
        matGen.fill_random_R (numCols, R_cur.data(), numRows, singularValues);
      }
    }

    /// \brief Generate a random test problem for the distributed-memory part of TSQR.
    ///
    /// Specifically, this routine generates an "R stack" test
    /// problem, where each (MPI) process has a square ncols by ncols
    /// upper triangular matrix A_local.  TSQR is supposed to factor
    /// the distributed matrix formed by "stacking" the A_local
    /// matrices on top of each other, so that MPI Rank 0's matrix is
    /// on top, Rank 1's matrix is below that, and so on.  (Ranks here
    /// are computed with respect to the given \c MessengerBase
    /// communicator wrapper.)
    ///
    /// \param A_local [out] ncols by ncols upper triangular matrix
    ///   (on each MPI process)
    ///
    /// \param A_global [out] Empty on all procs but Proc 0, where it
    ///   starts empty (not required, but this is more efficient) and
    ///   gets resized here
    ///
    /// \param ncols [in] Number of columns in the matrix to generate.
    ///   Number of rows in the matrix will be (# MPI processes)*ncols.
    ///
    /// \param generator [in/out] Normal(0,1) (pseudo)random number
    ///   generator
    ///
    /// \param messenger [in/out] MPI communicator object for Scalars
    ///
    template<class Ordinal, class Scalar, class Generator>
    void
    par_tsqr_test_problem (Generator& generator,
                           Matrix<Ordinal, Scalar>& A_local,
                           Matrix<Ordinal, Scalar>& A_global,
                           const Ordinal ncols,
                           const Teuchos::RCP<MessengerBase<Scalar> >& messenger)
    {
      const int nprocs = messenger->size();
      const int my_rank = messenger->rank();
      A_local.reshape (ncols, ncols);

      if (my_rank == 0)
        {
          typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;

          std::vector< magnitude_type > singular_values (ncols);
          singular_values[0] = magnitude_type(1);
          for (Ordinal k = 1; k < ncols; ++k)
            singular_values[k] = singular_values[k-1] / magnitude_type(2);

          generateStack (generator, A_global, &singular_values[0], nprocs, ncols);
          scatterStack (A_global, A_local, messenger);
        }
      else
        scatterStack (A_global, A_local, messenger);
    }


  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_generateStack_hpp
