// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_DualViewUtil.hpp"
#include "Teuchos_TestForException.hpp"

namespace Tpetra {
namespace Details {

auto view_alloc_no_init (const std::string& label) ->
  decltype (Kokkos::view_alloc (label, Kokkos::WithoutInitializing))
{
  return Kokkos::view_alloc (label, Kokkos::WithoutInitializing);
}

} // namespace Details
} // namespace Tpetra



