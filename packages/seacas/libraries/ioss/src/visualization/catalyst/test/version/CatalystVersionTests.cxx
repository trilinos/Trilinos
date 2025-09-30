// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Iovs_CatalystVersion.h>
#include <catch2/catch_test_macros.hpp>

TEST_CASE_METHOD(Iovs::CatalystVersion, "CatalystVersion", "[catalystVersion]")
{

  REQUIRE(getIOSSCatalystInterfaceVersion() == iossCatalystInterfaceVersion);
}

TEST_CASE_METHOD(Iovs::CatalystVersion, "PluginInterfaceSameVersion", "[catalystVersion]")
{

  REQUIRE(isInterfaceCompatibleWithPlugin("1.0.0", "1.0.0") == true);
  REQUIRE(isInterfaceCompatibleWithPlugin("1.3.0", "1.3.0") == true);
  REQUIRE(isInterfaceCompatibleWithPlugin("1.3.7", "1.3.7") == true);
}

TEST_CASE_METHOD(Iovs::CatalystVersion, "InterfaceOlderVersionThanPlugin", "[catalystVersion]")
{

  REQUIRE(isInterfaceCompatibleWithPlugin("1.0.0", "2.0.0") == false);
  REQUIRE(isInterfaceCompatibleWithPlugin("1.3.0", "1.4.0") == true);
  REQUIRE(isInterfaceCompatibleWithPlugin("1.5.72", "1.5.99") == true);
}

TEST_CASE_METHOD(Iovs::CatalystVersion, "InterfaceNewerVersionThanPlugin", "[catalystVersion]")
{

  REQUIRE(isInterfaceCompatibleWithPlugin("4.0.0", "1.0.0") == false);
  REQUIRE(isInterfaceCompatibleWithPlugin("2.7.2", "2.6.4") == false);
  REQUIRE(isInterfaceCompatibleWithPlugin("2.7.98", "2.7.90") == false);
}
