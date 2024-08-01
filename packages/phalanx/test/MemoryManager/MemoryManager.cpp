// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_MemoryManager.hpp"

#include "Sacado.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_View_Fad.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_DynRankView_Fad.hpp"

namespace phalanx_test {

  TEUCHOS_UNIT_TEST(Kokkos_AllocationTracker, MemoryManager)
  {
    PHX::MemoryManager pool1;

    auto pool2 = pool1.clone();
  }

}
