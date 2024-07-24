// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_CONFIGDEFS_HPP
#define TSQR_CONFIGDEFS_HPP

/// \file Tsqr_ConfigDefs.hpp
/// \brief Include this header to get TSQR's configuration options.

// Users should not include TpetraTSQR_config.h directly.
// Include Tsqr_ConfigDefs.hpp instead.
#include "TpetraTSQR_config.h"

/// \namespace TSQR
/// \brief Implementation of the Tall Skinny QR (TSQR) factorization.
///
/// This namespace contains a full hybrid-parallel (MPI + Kokkos)
/// implementation of the Tall Skinny QR (TSQR) factorization.  The
/// following paper describes the implementation:
///
/// Mark Hoemmen.  "A communication-avoiding, hybrid-parallel,
/// rank-revealing orthogonalization method."  IEEE International
/// Parallel and Distributed Processing Symposium (IPDPS), April 2011.
///
/// For further details, see the following:
///
/// Marghoob Mohiyuddin, Mark Hoemmen, James Demmel, and Kathy Yelick.
/// "Minimizing Communication in Sparse Matrix Solvers."  In
/// Proceedings of Supercomputing 2009, November 2009.
///
/// James Demmel, Laura Grigori, Mark Hoemmen, and Julien Langou.
/// "Communication-optimal parallel and sequential QR and LU
/// factorizations."  SIAM Journal on Scientific Computing, Volume 34,
/// Issue 1, 2012.
namespace TSQR {
  //
  // We declare the TSQR namespace here so that Doxygen will find it
  // and pull in all its documentation.
  //

  /// \namespace TSQR::Test
  /// \brief Accuracy and performance tests for TSQR.
  ///
  /// The classes and routines here are not intended for consumers of
  /// TSQR, but may be helpful as examples.
  namespace Test {
    //
    // We declare the TSQR::Test namespace here so that Doxygen will
    // find it and pull in all its documentation.
    //
  } // namespace Test
} // namespace TSQR

#endif // TSQR_CONFIGDEFS_HPP
