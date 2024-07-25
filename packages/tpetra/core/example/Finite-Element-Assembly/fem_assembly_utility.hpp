// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_UTILITY_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_UTILITY_HPP

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include <iostream>

namespace TpetraExamples {

template<typename V>
Kokkos::View<typename V::data_type, typename V::array_layout, typename V::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
makeUnmanaged(const V& v)
{
  return v;
}

// Return a pointer (RCP is like std::shared_ptr) to an output
// stream.  It prints on Process 0 of the given MPI communicator,
// but ignores all output on other MPI processes.
Teuchos::RCP<Teuchos::FancyOStream>
getOutputStream (const Teuchos::Comm<int>& comm)
{
  using Teuchos::getFancyOStream;

  const int myRank = comm.getRank ();
  if (0 == myRank) {
    // Process 0 of the given communicator prints to std::cout.
    return getFancyOStream (Teuchos::rcpFromRef (std::cout));
  }
  else {
    // A "black hole output stream" ignores all output directed to it.
    return getFancyOStream (Teuchos::rcp (new Teuchos::oblackholestream ()));
  }
}

} // namespace TpetraExamples

#endif  // TPETRAEXAMPLES_FEM_ASSEMBLY_UTILITY_HPP

