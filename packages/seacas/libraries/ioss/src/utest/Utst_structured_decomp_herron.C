// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Utst_structured_decomp.h"
#include <catch.hpp>

// Disable these tests on NVCC. It tries to optimize and takes forever to build...
#ifndef __NVCC__
TEST_CASE("herron-dutton", "[herron-dutton_zgc]")
{
  int                                       zone = 1;
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x128x728"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "A1", zones.back()->m_zone, "zone05", 5, Ioss::IJK_t{{3, 1, 2}}, Ioss::IJK_t{{129, 1, 1}},
      Ioss::IJK_t{{129, 129, 729}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{129, 729, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A2", zones.back()->m_zone, "zone04", 4, Ioss::IJK_t{{-3, -1, 2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 129, 729}}, Ioss::IJK_t{{129, 1, 1}}, Ioss::IJK_t{{1, 729, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A3", zones.back()->m_zone, "zone03", 3, Ioss::IJK_t{{2, 3, 1}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{129, 1, 729}}, Ioss::IJK_t{{1, 1, 193}}, Ioss::IJK_t{{729, 129, 193}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A4", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{-2, -3, 1}}, Ioss::IJK_t{{1, 129, 1}},
      Ioss::IJK_t{{129, 129, 729}}, Ioss::IJK_t{{1, 129, 193}}, Ioss::IJK_t{{729, 1, 193}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "728x128x192"));

  zones.back()->m_zoneConnectivity.emplace_back(
      "B4", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{3, -1, -2}}, Ioss::IJK_t{{1, 1, 193}},
      Ioss::IJK_t{{729, 129, 193}}, Ioss::IJK_t{{129, 129, 1}}, Ioss::IJK_t{{1, 129, 729}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A5", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{1, 3, -2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{729, 129, 1}}, Ioss::IJK_t{{129, 1, 1}}, Ioss::IJK_t{{857, 1, 129}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A6", zones.back()->m_zone, "zone05", 5, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{729, 1, 193}}, Ioss::IJK_t{{129, 1, 193}}, Ioss::IJK_t{{129, 729, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A7", zones.back()->m_zone, "zone04", 4, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{1, 129, 1}},
      Ioss::IJK_t{{729, 129, 193}}, Ioss::IJK_t{{1, 1, 193}}, Ioss::IJK_t{{1, 729, 1}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "728x128x192"));

  zones.back()->m_zoneConnectivity.emplace_back(
      "B3", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{3, 1, 2}}, Ioss::IJK_t{{1, 1, 193}},
      Ioss::IJK_t{{729, 129, 193}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{129, 1, 729}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A8", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{1, 3, -2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{729, 129, 1}}, Ioss::IJK_t{{129, 1, 1}}, Ioss::IJK_t{{857, 1, 129}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A9", zones.back()->m_zone, "zone04", 4, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{729, 1, 193}}, Ioss::IJK_t{{129, 1, 193}}, Ioss::IJK_t{{129, 729, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A10", zones.back()->m_zone, "zone05", 5, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{1, 129, 1}},
      Ioss::IJK_t{{729, 129, 193}}, Ioss::IJK_t{{1, 1, 193}}, Ioss::IJK_t{{1, 729, 1}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x728x192"));

  zones.back()->m_zoneConnectivity.emplace_back(
      "B2", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{-2, 3, -1}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{129, 729, 1}}, Ioss::IJK_t{{1, 129, 1}}, Ioss::IJK_t{{1, 1, 729}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "B7", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 729, 193}}, Ioss::IJK_t{{1, 129, 193}}, Ioss::IJK_t{{729, 129, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "B9", zones.back()->m_zone, "zone03", 3, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{129, 1, 1}},
      Ioss::IJK_t{{129, 729, 193}}, Ioss::IJK_t{{1, 1, 193}}, Ioss::IJK_t{{729, 1, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A11", zones.back()->m_zone, "zone07", 7, Ioss::IJK_t{{3, 1, 2}}, Ioss::IJK_t{{1, 1, 193}},
      Ioss::IJK_t{{129, 729, 193}}, Ioss::IJK_t{{129, 1, 1}}, Ioss::IJK_t{{857, 1, 129}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x728x192"));

  zones.back()->m_zoneConnectivity.emplace_back(
      "B1", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{2, 3, 1}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{129, 729, 1}}, Ioss::IJK_t{{129, 1, 1}}, Ioss::IJK_t{{129, 129, 729}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "B6", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{129, 1, 1}},
      Ioss::IJK_t{{129, 729, 193}}, Ioss::IJK_t{{1, 1, 193}}, Ioss::IJK_t{{729, 1, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "B10", zones.back()->m_zone, "zone03", 3, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 729, 193}}, Ioss::IJK_t{{1, 129, 193}}, Ioss::IJK_t{{729, 129, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A12", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{3, 1, 2}}, Ioss::IJK_t{{1, 1, 193}},
      Ioss::IJK_t{{129, 729, 193}}, Ioss::IJK_t{{129, 1, 1}}, Ioss::IJK_t{{857, 1, 129}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "856x192x128"));

  zones.back()->m_zoneConnectivity.emplace_back(
      "B8", zones.back()->m_zone, "zone03", 3, Ioss::IJK_t{{1, -3, 2}}, Ioss::IJK_t{{129, 1, 1}},
      Ioss::IJK_t{{857, 1, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{729, 129, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A13", zones.back()->m_zone, "zone07", 7, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{857, 193, 1}}, Ioss::IJK_t{{1, 1, 129}}, Ioss::IJK_t{{857, 193, 129}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A14", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 129}},
      Ioss::IJK_t{{857, 193, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{857, 193, 1}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "856x192x128"));

  zones.back()->m_zoneConnectivity.emplace_back(
      "B11", zones.back()->m_zone, "zone04", 4, Ioss::IJK_t{{2, 3, 1}}, Ioss::IJK_t{{129, 1, 1}},
      Ioss::IJK_t{{857, 1, 129}}, Ioss::IJK_t{{1, 1, 193}}, Ioss::IJK_t{{129, 729, 193}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "B13", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 129}},
      Ioss::IJK_t{{857, 193, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{857, 193, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A15", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{857, 193, 1}}, Ioss::IJK_t{{1, 1, 129}}, Ioss::IJK_t{{857, 193, 129}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "856x192x128"));

  zones.back()->m_zoneConnectivity.emplace_back(
      "B12", zones.back()->m_zone, "zone05", 5, Ioss::IJK_t{{2, 3, 1}}, Ioss::IJK_t{{129, 1, 1}},
      Ioss::IJK_t{{857, 1, 129}}, Ioss::IJK_t{{1, 1, 193}}, Ioss::IJK_t{{129, 729, 193}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "B14", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{857, 193, 1}}, Ioss::IJK_t{{1, 1, 129}}, Ioss::IJK_t{{857, 193, 129}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "A16", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 129}},
      Ioss::IJK_t{{857, 193, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{857, 193, 1}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "856x192x128"));

  zones.back()->m_zoneConnectivity.emplace_back(
      "B5", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{1, -3, 2}}, Ioss::IJK_t{{129, 1, 1}},
      Ioss::IJK_t{{857, 1, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{729, 129, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "B15", zones.back()->m_zone, "zone07", 7, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 129}},
      Ioss::IJK_t{{857, 193, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{857, 193, 1}});

  zones.back()->m_zoneConnectivity.emplace_back(
      "B16", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{857, 193, 1}}, Ioss::IJK_t{{1, 1, 129}}, Ioss::IJK_t{{857, 193, 129}});

  // There was an issue with smaller proc counts...
  double load_balance_tolerance = 1.4;
  for (size_t proc_count = 35; proc_count < 64; proc_count++) {
    std::string name = "Herron_Dutton_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count); }
  }

  load_balance_tolerance = 1.25;
  for (size_t proc_count = 4508; proc_count >= 10; proc_count /= 2) {
    std::string name = "Herron_Dutton_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count); }
  }

  cleanup(zones);
}
#endif
