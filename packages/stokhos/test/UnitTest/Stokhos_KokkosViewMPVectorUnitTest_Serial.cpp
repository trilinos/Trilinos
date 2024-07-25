// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_KokkosViewMPVectorUnitTest.hpp"

#include "Kokkos_Core.hpp"

// Instantiate test for Serial device
using Kokkos::Serial;
VIEW_MP_VECTOR_TESTS_DEVICE( Serial )
