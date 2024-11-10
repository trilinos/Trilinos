//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include "Tempus_Version.hpp"

namespace Tempus_Test {

TEUCHOS_UNIT_TEST(version, default) { Tempus::version(); }

}  // namespace Tempus_Test
