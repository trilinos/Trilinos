// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_RANDOM_DEF_HPP
#define TPETRA_DETAILS_RANDOM_DEF_HPP

#include "Tpetra_Details_Random_decl.hpp"
#include "Teuchos_TestForException.hpp"
#include "Kokkos_Random.hpp"

namespace Tpetra {
namespace Details {

namespace {  // anonymous
template <class ExecutionSpace>
void finalize_pool() {
  using PoolClass = Static_Random_XorShift64_Pool<ExecutionSpace>;
  if (PoolClass::pool_ != nullptr) {
    delete PoolClass::pool_;
  }
  PoolClass::pool_ = nullptr;
}
}  // end namespace

template <class ExecutionSpace>
unsigned int Static_Random_XorShift64_Pool<ExecutionSpace>::getSeedFromRank(int mpi_rank) {
  // Seed the pseudorandom number generator using the calling
  // process' rank.  This helps decorrelate different process'
  // pseudorandom streams.  It's not perfect but it's effective and
  // doesn't require MPI communication.  The seed also includes bits
  // from the standard library's rand().
  uint64_t myRank   = static_cast<uint64_t>(mpi_rank);
  uint64_t seed64   = static_cast<uint64_t>(std::rand()) + myRank + 17311uLL;
  unsigned int seed = static_cast<unsigned int>(seed64 & 0xffffffff);
  return seed;
}

template <class ExecutionSpace>
void Static_Random_XorShift64_Pool<ExecutionSpace>::
    resetPool(int mpi_rank) {
  if (isSet())
    delete pool_;
  else
    Kokkos::push_finalize_hook(finalize_pool<ExecutionSpace>);

  pool_ = new Kokkos::Random_XorShift64_Pool<ExecutionSpace>(getSeedFromRank(mpi_rank));
}

template <class ExecutionSpace>
bool Static_Random_XorShift64_Pool<ExecutionSpace>::
    isSet() {
  return pool_ != nullptr;
}

template <class ExecutionSpace>
Kokkos::Random_XorShift64_Pool<ExecutionSpace>& Static_Random_XorShift64_Pool<ExecutionSpace>::
    getPool() {
  TEUCHOS_TEST_FOR_EXCEPTION(!isSet(), std::runtime_error, "Tpetra::Details::Static_Random_XorShift64_Pool: resetPool() must be called before getPool");
  return *pool_;
}

}  // namespace Details
}  // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
// NOTE: Tpetra::Details::Random is templated on execution space,
// but we can use this trick to get it out of the node.
//

#define TPETRA_DETAILS_RANDOM_INSTANT(NODE)                                                                                                                \
  template <>                                                                                                                                              \
  Kokkos::Random_XorShift64_Pool<typename NODE::execution_space>* Details::Static_Random_XorShift64_Pool<typename NODE::execution_space>::pool_ = nullptr; \
  template class Details::Static_Random_XorShift64_Pool<typename NODE::execution_space>;

#endif  // TPETRA_DETAILS_RANDOM_DEF_HPP
