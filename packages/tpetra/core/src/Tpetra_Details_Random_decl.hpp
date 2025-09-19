// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_RANDOM_DECL_HPP
#define TPETRA_DETAILS_RANDOM_DECL_HPP

#include "TpetraCore_config.h"
#include "Kokkos_Random.hpp"

namespace Tpetra {
namespace Details {

template <class ExecutionSpace>
class Static_Random_XorShift64_Pool {
 public:
  // The resetPool function will re-initialize the pool based on the system RNG and the MPI rank.
  // On GPU architectures, this will likely involve non-trivial host-to-device transfers.
  static void resetPool(int mpi_rank);

  // The isSet function returns true if resetPool has been callled.
  static bool isSet();
  // The getPool function will return the existing pool.
  static Kokkos::Random_XorShift64_Pool<ExecutionSpace>& getPool();

  // Do not access this directly.  This is public only for deallocation purposes
  static Kokkos::Random_XorShift64_Pool<ExecutionSpace>* pool_;

 private:
  static unsigned int getSeedFromRank(int mpi_rank);
};

}  // namespace Details
}  // namespace Tpetra

#endif  // TPETRA_DETAILS_RANDOM_DECL_HPP
