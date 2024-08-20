// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_ZoneConnectivity.h"
#include "cgns/Iocgns_StructuredZoneData.h"
#include "cgns/Iocgns_Utils.h"
#include <assert.h>
#include <catch2/catch_test_macros.hpp>
#include <map>
#include <numeric>
#include <set>
#include <stddef.h>
#include <stdint.h>
#include <string>
#include <vector>

#include "Ioss_CodeTypes.h"
#include "Ioss_Utils.h"

namespace {
  int64_t generate_guid(size_t id, int rank, int proc_count)
  {
    static size_t lpow2 = 0;
    static int    procs = -1;
    if (procs != proc_count) {
      lpow2 = Ioss::Utils::log_power_2(proc_count);
      procs = proc_count;
    }
    assert(rank >= 0);
    return (id << lpow2) + rank;
  }

  void update_zgc_data(std::vector<Iocgns::StructuredZoneData *> &zones, int proc_count)
  {
    for (auto &zone : zones) {
      if (zone->is_active()) {
        zone->resolve_zgc_split_donor(zones);
      }
    }

    for (auto &zone : zones) {
      if (zone->is_active()) {
        zone->update_zgc_processor(zones);
      }
    }

    // Set GUID on all zgc instances...
    for (auto &zone : zones) {
      if (zone->is_active()) {
        for (auto &zgc : zone->m_zoneConnectivity) {
          zgc.m_ownerGUID = generate_guid(zgc.m_ownerZone, zgc.m_ownerProcessor, proc_count);
          zgc.m_donorGUID = generate_guid(zgc.m_donorZone, zgc.m_donorProcessor, proc_count);
        }
      }
    }
  }

} // namespace

void cleanup(std::vector<Iocgns::StructuredZoneData *> &zones)
{
  for (auto &zone : zones) {
    delete zone;
    zone = nullptr;
  }
}

// Helper function that drives all tests...
void check_split_assign(std::vector<Iocgns::StructuredZoneData *> &zones,
                        double load_balance_tolerance, size_t proc_count, double min_toler = 0.9,
                        double max_toler = 1.0)
{
#if IOSS_DEBUG_OUTPUT
  bool verbose = true;
#else
  bool verbose = false;
#endif

  double total_work =
      std::accumulate(zones.begin(), zones.end(), 0.0,
                      [](double a, Iocgns::StructuredZoneData *b) { return a + b->work(); });

  double avg_work = total_work / (double)proc_count;
  SECTION("split_zones")
  {
    Iocgns::Utils::pre_split(zones, avg_work, load_balance_tolerance, 0, proc_count, verbose);

    double max_work = avg_work * load_balance_tolerance * max_toler;
    for (const auto zone : zones) {
      if (zone->is_active()) {
        CHECK(zone->work() <= max_work);
      }
    }

    for (size_t i = 0; i < zones.size(); i++) {
      CHECK(zones[i]->m_zone == int(i) + 1);
    }

    SECTION("assign_to_procs")
    {
      std::vector<size_t> work_vector(proc_count);
      Iocgns::Utils::assign_zones_to_procs(zones, work_vector, verbose);

#if 0
        fmt::print(stderr, "\nDecomposition for {} processors; Total work = {}, Average = {}\n",
                   proc_count, fmt::group_digits((size_t)total_work), fmt::group_digits((size_t)avg_work));

          for (const auto zone : zones) {
            if (zone->is_active()) {
              fmt::print(stderr, "Zone {}\tProc: {}\tOrdinal: {}x{}x{}\tWork: {}\n",
                         zone->m_name, zone->m_proc, zone->m_ordinal[0], zone->m_ordinal[1],
                         zone->m_ordinal[2], fmt::group_digits(zone->work()));
            }
          }
#endif
      // Each active zone must be on a processor
      for (const auto zone : zones) {
        if (zone->is_active()) {
          CHECK(zone->m_proc >= 0);
        }
      }

      // Work must be min_work <= work <= max_work
      double min_work = avg_work / load_balance_tolerance * min_toler;
      for (auto work : work_vector) {
        CHECK(work >= min_work);
        CHECK(work <= max_work * max_toler);
      }

      // A processor cannot have more than one zone with the same adam zone
      std::set<std::pair<int, int>> proc_adam_map;
      for (const auto zone : zones) {
        if (zone->is_active()) {
          auto success = proc_adam_map.insert(std::make_pair(zone->m_adam->m_zone, zone->m_proc));
          CHECK(success.second);
        }
      }

      // Zone Grid Connectivity Checks:
      update_zgc_data(zones, proc_count);

      // Zone Grid Connectivity instances can't connect to themselves...
      for (auto &zone : zones) {
        if (zone->is_active()) {
          for (const auto &zgc : zone->m_zoneConnectivity) {
            if (zgc.is_active()) {
              CHECK(zgc.m_ownerZone != zgc.m_donorZone);
              CHECK(zgc.m_ownerGUID != zgc.m_donorGUID);
            }
          }
        }
      }

      // In Iocgns::Utils::common_write_metadata, there is code to make
      // sure that the zgc.m_connectionName  is unique for all zgc instances on
      // a zone / processor pair (if !parallel_io which is file-per-processor)
      // The uniquification appends a letter from 'A' to 'Z' to the name
      // If the name is still not unique, repeats process with 'AA' to 'ZZ'
      // Make sure that there are not more than 26 + 26*26 + 1 instances of the same
      // name on a zone to ensure that this works...
      for (auto &zone : zones) {
        if (zone->is_active()) {
          std::map<std::string, int> zgc_map;
          for (const auto &zgc : zone->m_zoneConnectivity) {
            if (zgc.is_active() && !zgc.is_from_decomp()) {
              zgc_map[zgc.m_connectionName]++;
            }
          }
          for (const auto &kk : zgc_map) {
            CHECK(kk.second < 26 * 26 + 26 + 1);
          }
        }
      } //

      // Zone Grid Connectivity from_decomp instances must be symmetric...
      // The GUID encodes the id and the processor,
      std::map<std::pair<size_t, size_t>, int> is_symm;
      for (auto &zone : zones) {
        if (zone->is_active()) {
          for (const auto &zgc : zone->m_zoneConnectivity) {
            if (zgc.is_active()) {
              is_symm[std::make_pair(std::min(zgc.m_ownerGUID, zgc.m_donorGUID),
                                     std::max(zgc.m_ownerGUID, zgc.m_donorGUID))]++;
            }
          }
        }
      }
      // Iterate `is_symm` and make sure there is an even number for all entries.
      for (const auto &item : is_symm) {
        CHECK(item.second % 2 == 0);
      }
    }
  }
}

TEST_CASE("single block")
{
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(1, "4x4x1"));

  int    proc_count             = 2;
  double load_balance_tolerance = 1.2;

  check_split_assign(zones, load_balance_tolerance, proc_count);
  cleanup(zones);
}

TEST_CASE("single block line")
{
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(1, "40x40x1"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::I;

  int    proc_count             = 4;
  double load_balance_tolerance = 1.05;

  check_split_assign(zones, load_balance_tolerance, proc_count);
  cleanup(zones);
}

TEST_CASE("single block ij-line")
{
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(1, "40x40x40"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::I | Iocgns::Ordinal::J;

  int    proc_count             = 4;
  double load_balance_tolerance = 1.05;

  check_split_assign(zones, load_balance_tolerance, proc_count);
  cleanup(zones);
}

TEST_CASE("single block jk-line")
{
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(1, "40x40x40"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J | Iocgns::Ordinal::K;

  int    proc_count             = 4;
  double load_balance_tolerance = 1.05;

  check_split_assign(zones, load_balance_tolerance, proc_count);
  cleanup(zones);
}

TEST_CASE("single block ik-line")
{
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(1, "40x40x40"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::I | Iocgns::Ordinal::K;

  int    proc_count             = 4;
  double load_balance_tolerance = 1.05;

  check_split_assign(zones, load_balance_tolerance, proc_count);
  cleanup(zones);
}

TEST_CASE("cube_2blocks")
{
  int                                       zone = 1;
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "50x20x50"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "A1", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{1, -3, 2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{6, 1, 6}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{6, 6, 1}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "50x50x30"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B1", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{1, 3, -2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{6, 6, 1}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{6, 1, 6}});
  double load_balance_tolerance = 1.2;

  for (size_t proc_count = 2; proc_count < 8; proc_count += 2) {
    std::string name = "cube_2blocks_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count, 0.9, 1.1); }
  }
  cleanup(zones);
}

TEST_CASE("2blocks_example")
{
  int                                       zone = 1;
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x4x1"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x8x1"));
  double load_balance_tolerance = 1.1;

  check_split_assign(zones, load_balance_tolerance, 4, 0.9, 1.1);
  cleanup(zones);
}

TEST_CASE("bump")
{
  int                                       zone = 1;
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "2x2x1"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "A1", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{3, 1, 1}},
      Ioss::IJK_t{{3, 3, 2}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{1, 3, 2}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "2x2x1"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B1", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 3, 2}}, Ioss::IJK_t{{3, 1, 1}}, Ioss::IJK_t{{3, 3, 2}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A2", zones.back()->m_zone, "zone03", 3, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{3, 1, 1}},
      Ioss::IJK_t{{3, 3, 2}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{1, 3, 2}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "2x2x1"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B2", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 3, 2}}, Ioss::IJK_t{{3, 1, 1}}, Ioss::IJK_t{{3, 3, 2}});
  double load_balance_tolerance = 1.2;

  SECTION("bump_ProcCount_2") { check_split_assign(zones, load_balance_tolerance, 2, 0.8, 1.2); }
  cleanup(zones);
}

TEST_CASE("bump_loose")
{
  int                                       zone = 1;
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "2x2x1"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "A1", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{3, 1, 1}},
      Ioss::IJK_t{{3, 3, 2}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{1, 3, 2}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "2x2x1"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B1", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 3, 2}}, Ioss::IJK_t{{3, 1, 1}}, Ioss::IJK_t{{3, 3, 2}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A2", zones.back()->m_zone, "zone03", 3, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{3, 1, 1}},
      Ioss::IJK_t{{3, 3, 2}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{1, 3, 2}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "2x2x1"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B2", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 3, 2}}, Ioss::IJK_t{{3, 1, 1}}, Ioss::IJK_t{{3, 3, 2}});
  double load_balance_tolerance = 1.4;

  SECTION("bump_loose_ProcCount_2")
  {
    check_split_assign(zones, load_balance_tolerance, 2, 0.8, 1.2);
  }
  cleanup(zones);
}

TEST_CASE("farmer plenum")
{
  std::vector<Iocgns::StructuredZoneData *> zones;
  int                                       zone = 1;

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "56x128x48"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "A1", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{2, -1, 3}}, Ioss::IJK_t{{57, 1, 1}},
      Ioss::IJK_t{{57, 33, 49}}, Ioss::IJK_t{{33, 1, 1}}, Ioss::IJK_t{{1, 1, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A2", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{57, 33, 1}},
      Ioss::IJK_t{{57, 65, 49}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{1, 33, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A3", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{57, 65, 1}},
      Ioss::IJK_t{{57, 97, 49}}, Ioss::IJK_t{{1, 33, 1}}, Ioss::IJK_t{{1, 65, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A4", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{-2, 1, 3}}, Ioss::IJK_t{{57, 97, 1}},
      Ioss::IJK_t{{57, 129, 49}}, Ioss::IJK_t{{1, 65, 1}}, Ioss::IJK_t{{33, 65, 49}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x64x48"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B1", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{-2, 1, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{33, 1, 49}}, Ioss::IJK_t{{57, 33, 1}}, Ioss::IJK_t{{57, 1, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B2", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 33, 49}}, Ioss::IJK_t{{57, 33, 1}}, Ioss::IJK_t{{57, 65, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B3", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 33, 1}},
      Ioss::IJK_t{{1, 65, 49}}, Ioss::IJK_t{{57, 65, 1}}, Ioss::IJK_t{{57, 97, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B4", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{2, -1, 3}}, Ioss::IJK_t{{1, 65, 1}},
      Ioss::IJK_t{{33, 65, 49}}, Ioss::IJK_t{{57, 97, 1}}, Ioss::IJK_t{{57, 129, 49}});

  double load_balance_tolerance = 1.1;

  for (size_t proc_count = 2; proc_count < 20; proc_count++) {
    std::string name = "Plenum_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count, 0.75); }
  }
  cleanup(zones);
}

TEST_CASE("grv-nose")
{
  std::vector<Iocgns::StructuredZoneData *> zones;
  int                                       zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "16x16x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "A1", zones.back()->m_zone, "zone4", 5, Ioss::IJK_t{{-2, -1, -3}}, Ioss::IJK_t{{1, 17, 1}},
      Ioss::IJK_t{{17, 17, 129}}, Ioss::IJK_t{{33, 17, 129}}, Ioss::IJK_t{{33, 1, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A2", zones.back()->m_zone, "zone2", 3, Ioss::IJK_t{{2, -1, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{17, 1, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{1, 17, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A3", zones.back()->m_zone, "zone5", 6, Ioss::IJK_t{{-1, 2, -3}}, Ioss::IJK_t{{17, 1, 1}},
      Ioss::IJK_t{{17, 17, 129}}, Ioss::IJK_t{{33, 1, 129}}, Ioss::IJK_t{{33, 17, 1}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x16x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B11", zones.back()->m_zone, "zone5", 6, Ioss::IJK_t{{-1, -2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{1, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B14", zones.back()->m_zone, "zone8", 9, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{65, 1, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{65, 17, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B15", zones.back()->m_zone, "zone9", 10, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 17, 1}},
      Ioss::IJK_t{{65, 17, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{65, 1, 129}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "16x16x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B2", zones.back()->m_zone, "zone1", 1, Ioss::IJK_t{{-2, 1, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{17, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A4", zones.back()->m_zone, "zone3", 4, Ioss::IJK_t{{-2, -1, -3}}, Ioss::IJK_t{{1, 17, 1}},
      Ioss::IJK_t{{17, 17, 129}}, Ioss::IJK_t{{33, 17, 129}}, Ioss::IJK_t{{33, 1, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A5", zones.back()->m_zone, "zone6", 7, Ioss::IJK_t{{1, -2, -3}}, Ioss::IJK_t{{17, 1, 1}},
      Ioss::IJK_t{{17, 17, 129}}, Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{1, 1, 1}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x16x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B4", zones.back()->m_zone, "zone2", 3, Ioss::IJK_t{{-2, -1, -3}}, Ioss::IJK_t{{33, 1, 1}},
      Ioss::IJK_t{{33, 17, 129}}, Ioss::IJK_t{{17, 17, 129}}, Ioss::IJK_t{{1, 17, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A6", zones.back()->m_zone, "zone6", 7, Ioss::IJK_t{{-1, -2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{33, 1, 129}}, Ioss::IJK_t{{33, 1, 1}}, Ioss::IJK_t{{1, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A7", zones.back()->m_zone, "zone9", 10, Ioss::IJK_t{{-1, -2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{1, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A8", zones.back()->m_zone, "zone5", 6, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 17, 1}},
      Ioss::IJK_t{{33, 17, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{33, 1, 129}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x16x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B1", zones.back()->m_zone, "zone1", 1, Ioss::IJK_t{{-2, -1, -3}}, Ioss::IJK_t{{33, 1, 1}},
      Ioss::IJK_t{{33, 17, 129}}, Ioss::IJK_t{{17, 17, 129}}, Ioss::IJK_t{{1, 17, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A9", zones.back()->m_zone, "zone5", 6, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{33, 1, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{33, 17, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A10", zones.back()->m_zone, "zone8", 9, Ioss::IJK_t{{-1, -2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{1, 1, 129}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x16x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B3", zones.back()->m_zone, "zone1", 1, Ioss::IJK_t{{-1, 2, -3}}, Ioss::IJK_t{{33, 1, 1}},
      Ioss::IJK_t{{33, 17, 129}}, Ioss::IJK_t{{17, 1, 129}}, Ioss::IJK_t{{17, 17, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B8", zones.back()->m_zone, "zone3", 4, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{33, 1, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{33, 17, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B9", zones.back()->m_zone, "zone4", 5, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 17, 1}},
      Ioss::IJK_t{{33, 17, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{33, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A11", zones.back()->m_zone, "zone10", 2, Ioss::IJK_t{{-1, -2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{1, 1, 129}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x16x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B5", zones.back()->m_zone, "zone2", 3, Ioss::IJK_t{{1, -2, -3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{17, 17, 129}}, Ioss::IJK_t{{17, 1, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B6", zones.back()->m_zone, "zone3", 4, Ioss::IJK_t{{-1, -2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{33, 1, 129}}, Ioss::IJK_t{{33, 1, 1}}, Ioss::IJK_t{{1, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A12", zones.back()->m_zone, "zone7", 8, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{33, 1, 1}},
      Ioss::IJK_t{{33, 17, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{1, 17, 129}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x16x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B12", zones.back()->m_zone, "zone6", 7, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{33, 1, 1}}, Ioss::IJK_t{{33, 17, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A13", zones.back()->m_zone, "zone9", 10, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{65, 1, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{65, 17, 129}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x16x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B10", zones.back()->m_zone, "zone4", 5, Ioss::IJK_t{{-1, -2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{1, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A14", zones.back()->m_zone, "zone10", 2, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 17, 1}},
      Ioss::IJK_t{{65, 17, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{65, 1, 129}});
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x16x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B7", zones.back()->m_zone, "zone3", 4, Ioss::IJK_t{{-1, -2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{1, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B13", zones.back()->m_zone, "zone7", 8, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 17, 1}},
      Ioss::IJK_t{{65, 17, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{65, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A15", zones.back()->m_zone, "zone10", 2, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{65, 1, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{65, 17, 129}});

  for (size_t proc_count = 3; proc_count <= 384; proc_count *= 2) {
    std::string name = "GRV-Nose_ProcCount_" + std::to_string(proc_count);
    SECTION(name)
    {
      double load_balance_tolerance = 1.2;
      check_split_assign(zones, load_balance_tolerance, proc_count, 0.9, 1.2);
    }
  }
  cleanup(zones);
}

TEST_CASE("grv")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x1x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x1x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x4x4"));

  for (size_t proc_count = 2; proc_count < 16; proc_count++) {
    std::string name = "GRV_ProcCount_" + std::to_string(proc_count);
    SECTION(name)
    {
      double load_balance_tolerance = 1.3;
      check_split_assign(zones, load_balance_tolerance, proc_count, .7, 1.1);
    }
  }
  cleanup(zones);
}

TEST_CASE("grv-large")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x16x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x16x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));

  for (size_t proc_count = 2; proc_count < 8192; proc_count *= 2) {
    std::string name = "GRV-LARGE_ProcCount_" + std::to_string(proc_count);
    SECTION(name)
    {
      double load_balance_tolerance = 1.3;
      check_split_assign(zones, load_balance_tolerance, proc_count, .7);
    }
  }
  cleanup(zones);
}

TEST_CASE("grv-large-ordinal")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x16x64"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x32x32"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x16x64"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x64x64"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;

  for (size_t proc_count = 2; proc_count < 8192; proc_count *= 2) {
    std::string name = "GRV-LARGE_ORDINAL_ProcCount_" + std::to_string(proc_count);
    SECTION(name)
    {
      double load_balance_tolerance = 1.3;
      check_split_assign(zones, load_balance_tolerance, proc_count, .7);
    }
  }
  cleanup(zones);
}

TEST_CASE("mk21")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x2x8"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x2x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x2x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x4x4"));

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x8x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "6x8x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x4"));

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x4x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x4x2"));

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4x2x2"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "16x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "16x4x4"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "16x4x4"));

  double load_balance_tolerance = 1.2;

  for (size_t proc_count = 2; proc_count < 17; proc_count++) {
    std::string name = "MK21_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count); }
  }
  cleanup(zones);
}

TEST_CASE("mk21-large")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x32x128"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x128x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x128x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x64x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x64x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "256x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "256x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "256x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x32x128"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x128x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "96x128x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x64x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x64x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "64x32x32"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "256x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "256x64x64"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "256x64x64"));

  double load_balance_tolerance = 1.1;

  for (size_t proc_count = 2; proc_count < 257; proc_count *= 2) {
    std::string name = "MK21_Large_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count, .8, 1.2); }
  }
  cleanup(zones);
}

TEST_CASE("farmer_h1_nozzle")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  // tried a nozzle only case:
  // farmer_h1_nozzle_3dh_2M-hex-join-domains.cgns. This
  // case runs on 320 processes, but on 384 processes, it results in
  // asymmetric communication, which appears to be coming from
  // IOSS.
  //
  // For example, on process 067, there are two interfaces
  // connected to 309:
  // current surf 189 1to1ConnectionB7 ranks: 67 309 gids: 2115 1845 owner ranges: 1 1 2 33 80 81
  // current surf 209 1to1ConnectionB7 ranks: 67 309 gids: 2115 1845 owner ranges: 1 1 2 33 70 71
  //
  // But on process 309, there is only one surface connection to 67:
  // current surf 4 1to1ConnectionA7 ranks: 309 67 gids: 1845 2115 owner ranges: 57 57 59 60 34 65

  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "56x128x48"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x64x48"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "56x192x128"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x64x128"));

  double load_balance_tolerance = 1.33;

  for (size_t proc_count = 3; proc_count <= 384; proc_count *= 2) {
    std::string name = "NOZ_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count); }
  }
  cleanup(zones);
}

TEST_CASE("farmer_h1_mk21")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  // We are still having some issues with IOSS auto decomp on the full
  // grid.
  // When trying to run it on 320 processes, I get the error that the
  // block global ids of an active interface are equivalent on both
  // sides of the interface. Combing through the pointwise we didn’t
  // find any connection that has two blocks the same, so it seems to
  // come from decomposition.
  // The specific connection that fails 6_293--6_294 on proc 316. Global ids are 3388 3388.

  // StructuredBlock 'blk-01' 192x64x88     1081344 cells,      1116505 nodes
  // StructuredBlock 'blk-02' 32x88x64      180224 cells,       190905 nodes
  // StructuredBlock 'blk-03' 32x64x88      180224 cells,       190905 nodes
  // StructuredBlock 'blk-04' 56x128x48      344064 cells,       360297 nodes
  // StructuredBlock 'blk-05' 32x64x48       98304 cells,       105105 nodes
  // StructuredBlock 'blk-06' 56x192x128     1376256 cells,      1419129 nodes
  // StructuredBlock 'blk-07' 72x128x64      589824 cells,       612105 nodes
  // StructuredBlock 'blk-08' 56x80x32      143360 cells,       152361 nodes
  // StructuredBlock 'blk-09' 8x48x128       49152 cells,        56889 nodes
  // StructuredBlock 'blk-10' 128x16x16       32768 cells,        37281 nodes
  // StructuredBlock 'blk-11' 32x64x128      262144 cells,       276705 nodes
  // StructuredBlock 'blk-12' 56x8x64       28672 cells,        33345 nodes

  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "192x64x88"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "A1", zones.back()->m_zone, "zone3", 3, Ioss::IJK_t{{2, -1, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{65, 1, 89}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{1, 65, 89}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A2", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{-3, -1, 2}}, Ioss::IJK_t{{129, 1, 1}},
      Ioss::IJK_t{{193, 1, 89}}, Ioss::IJK_t{{1, 1, 65}}, Ioss::IJK_t{{1, 89, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A3", zones.back()->m_zone, "zone07", 7, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{65, 1, 1}},
      Ioss::IJK_t{{129, 1, 65}}, Ioss::IJK_t{{73, 33, 65}}, Ioss::IJK_t{{73, 97, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A4", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{-3, 2, 1}}, Ioss::IJK_t{{65, 1, 65}},
      Ioss::IJK_t{{129, 1, 73}}, Ioss::IJK_t{{1, 49, 97}}, Ioss::IJK_t{{9, 49, 33}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A5", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{3, 1, 2}}, Ioss::IJK_t{{65, 1, 73}},
      Ioss::IJK_t{{129, 1, 89}}, Ioss::IJK_t{{57, 177, 33}}, Ioss::IJK_t{{57, 193, 97}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x88x64"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B2", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{-2, 3, -1}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 89, 65}}, Ioss::IJK_t{{193, 1, 1}}, Ioss::IJK_t{{129, 1, 89}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A6", zones.back()->m_zone, "zone07", 7, Ioss::IJK_t{{2, -3, -1}}, Ioss::IJK_t{{1, 1, 65}},
      Ioss::IJK_t{{33, 65, 65}}, Ioss::IJK_t{{73, 97, 65}}, Ioss::IJK_t{{73, 129, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A7", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{-3, 1, -2}}, Ioss::IJK_t{{1, 65, 65}},
      Ioss::IJK_t{{33, 73, 65}}, Ioss::IJK_t{{1, 49, 33}}, Ioss::IJK_t{{9, 49, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A8", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{3, 2, -1}}, Ioss::IJK_t{{1, 73, 65}},
      Ioss::IJK_t{{33, 89, 65}}, Ioss::IJK_t{{57, 177, 97}}, Ioss::IJK_t{{57, 193, 129}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x64x88"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B1", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{-2, 1, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 65, 89}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{65, 1, 89}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A9", zones.back()->m_zone, "zone07", 7, Ioss::IJK_t{{-2, -1, -3}}, Ioss::IJK_t{{1, 65, 1}},
      Ioss::IJK_t{{33, 65, 65}}, Ioss::IJK_t{{73, 33, 65}}, Ioss::IJK_t{{73, 1, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A10", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{3, -2, 1}}, Ioss::IJK_t{{1, 65, 65}},
      Ioss::IJK_t{{33, 65, 73}}, Ioss::IJK_t{{1, 49, 97}}, Ioss::IJK_t{{9, 49, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A11", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{-3, -1, 2}}, Ioss::IJK_t{{1, 65, 73}},
      Ioss::IJK_t{{33, 65, 89}}, Ioss::IJK_t{{57, 177, 33}}, Ioss::IJK_t{{57, 193, 1}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "56x128x48"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "A12", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{1, 3, -2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{57, 129, 1}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{57, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A13", zones.back()->m_zone, "zone05", 5, Ioss::IJK_t{{2, -1, 3}}, Ioss::IJK_t{{57, 1, 1}},
      Ioss::IJK_t{{57, 33, 49}}, Ioss::IJK_t{{33, 1, 1}}, Ioss::IJK_t{{1, 1, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A14", zones.back()->m_zone, "zone05", 5, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{57, 33, 1}},
      Ioss::IJK_t{{57, 97, 49}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{1, 65, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A15", zones.back()->m_zone, "zone05", 5, Ioss::IJK_t{{-2, 1, 3}}, Ioss::IJK_t{{57, 97, 1}},
      Ioss::IJK_t{{57, 129, 49}}, Ioss::IJK_t{{1, 65, 1}}, Ioss::IJK_t{{33, 65, 49}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x64x48"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B13", zones.back()->m_zone, "zone04", 4, Ioss::IJK_t{{-2, 1, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{33, 1, 49}}, Ioss::IJK_t{{57, 33, 1}}, Ioss::IJK_t{{57, 1, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B14", zones.back()->m_zone, "zone04", 4, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 65, 49}}, Ioss::IJK_t{{57, 33, 1}}, Ioss::IJK_t{{57, 97, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B15", zones.back()->m_zone, "zone04", 4, Ioss::IJK_t{{2, -1, 3}}, Ioss::IJK_t{{1, 65, 1}},
      Ioss::IJK_t{{33, 65, 49}}, Ioss::IJK_t{{57, 97, 1}}, Ioss::IJK_t{{57, 129, 49}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A16", zones.back()->m_zone, "zone11", 11, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{33, 65, 1}}, Ioss::IJK_t{{1, 1, 129}}, Ioss::IJK_t{{33, 65, 129}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "56x192x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B5", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{2, 3, 1}}, Ioss::IJK_t{{57, 177, 33}},
      Ioss::IJK_t{{57, 193, 97}}, Ioss::IJK_t{{65, 1, 73}}, Ioss::IJK_t{{129, 1, 89}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B8", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{-3, 2, 1}}, Ioss::IJK_t{{57, 177, 97}},
      Ioss::IJK_t{{57, 193, 129}}, Ioss::IJK_t{{1, 73, 65}}, Ioss::IJK_t{{33, 89, 65}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B11", zones.back()->m_zone, "zone03", 3, Ioss::IJK_t{{-2, 3, -1}}, Ioss::IJK_t{{57, 177, 1}},
      Ioss::IJK_t{{57, 193, 33}}, Ioss::IJK_t{{33, 65, 73}}, Ioss::IJK_t{{1, 65, 89}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B12", zones.back()->m_zone, "zone04", 4, Ioss::IJK_t{{1, -3, 2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{57, 1, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{57, 129, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A17", zones.back()->m_zone, "zone11", 11, Ioss::IJK_t{{-2, -3, 1}}, Ioss::IJK_t{{57, 1, 97}},
      Ioss::IJK_t{{57, 129, 129}}, Ioss::IJK_t{{1, 65, 129}}, Ioss::IJK_t{{33, 65, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A18", zones.back()->m_zone, "zone11", 11, Ioss::IJK_t{{1, -3, 2}}, Ioss::IJK_t{{57, 1, 33}},
      Ioss::IJK_t{{57, 129, 97}}, Ioss::IJK_t{{1, 1, 129}}, Ioss::IJK_t{{1, 65, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A19", zones.back()->m_zone, "zone11", 11, Ioss::IJK_t{{2, -3, -1}}, Ioss::IJK_t{{57, 1, 1}},
      Ioss::IJK_t{{57, 129, 33}}, Ioss::IJK_t{{33, 1, 129}}, Ioss::IJK_t{{1, 1, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A20", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{-1, 2, -3}}, Ioss::IJK_t{{57, 129, 1}},
      Ioss::IJK_t{{57, 177, 129}}, Ioss::IJK_t{{9, 1, 129}}, Ioss::IJK_t{{9, 49, 1}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "72x128x64"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B3", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{73, 33, 1}},
      Ioss::IJK_t{{73, 97, 65}}, Ioss::IJK_t{{65, 1, 65}}, Ioss::IJK_t{{129, 1, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B6", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{-3, 1, -2}}, Ioss::IJK_t{{73, 97, 1}},
      Ioss::IJK_t{{73, 129, 65}}, Ioss::IJK_t{{1, 65, 65}}, Ioss::IJK_t{{33, 1, 65}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B9", zones.back()->m_zone, "zone03", 3, Ioss::IJK_t{{-2, -1, -3}}, Ioss::IJK_t{{73, 1, 1}},
      Ioss::IJK_t{{73, 33, 65}}, Ioss::IJK_t{{33, 65, 65}}, Ioss::IJK_t{{1, 65, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A21", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{1, 3, -2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{57, 33, 1}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{57, 1, 33}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A22", zones.back()->m_zone, "zone10", 10, Ioss::IJK_t{{2, -1, 3}}, Ioss::IJK_t{{57, 1, 1}},
      Ioss::IJK_t{{73, 129, 1}}, Ioss::IJK_t{{129, 1, 17}}, Ioss::IJK_t{{1, 17, 17}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A23", zones.back()->m_zone, "zone12", 12, Ioss::IJK_t{{1, 3, -2}}, Ioss::IJK_t{{1, 33, 1}},
      Ioss::IJK_t{{57, 97, 1}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{57, 1, 65}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A24", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{1, -3, 2}}, Ioss::IJK_t{{1, 97, 1}},
      Ioss::IJK_t{{57, 129, 1}}, Ioss::IJK_t{{1, 81, 33}}, Ioss::IJK_t{{57, 81, 1}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "56x80x32"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B21", zones.back()->m_zone, "zone07", 7, Ioss::IJK_t{{1, -3, 2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{57, 1, 33}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{57, 33, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B24", zones.back()->m_zone, "zone07", 7, Ioss::IJK_t{{1, 3, -2}}, Ioss::IJK_t{{1, 81, 1}},
      Ioss::IJK_t{{57, 81, 33}}, Ioss::IJK_t{{1, 129, 1}}, Ioss::IJK_t{{57, 97, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A25", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{57, 1, 1}},
      Ioss::IJK_t{{57, 9, 33}}, Ioss::IJK_t{{1, 1, 129}}, Ioss::IJK_t{{9, 1, 97}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A26", zones.back()->m_zone, "zone11", 11, Ioss::IJK_t{{3, 2, -1}}, Ioss::IJK_t{{57, 9, 1}},
      Ioss::IJK_t{{57, 73, 33}}, Ioss::IJK_t{{33, 1, 1}}, Ioss::IJK_t{{1, 65, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A27", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{2, -1, 3}}, Ioss::IJK_t{{57, 73, 1}},
      Ioss::IJK_t{{57, 81, 33}}, Ioss::IJK_t{{9, 1, 1}}, Ioss::IJK_t{{1, 1, 33}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A28", zones.back()->m_zone, "zone12", 12, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 33}},
      Ioss::IJK_t{{57, 9, 33}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{57, 9, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A29", zones.back()->m_zone, "zone12", 12, Ioss::IJK_t{{1, 3, -2}}, Ioss::IJK_t{{1, 9, 33}},
      Ioss::IJK_t{{57, 73, 33}}, Ioss::IJK_t{{1, 9, 1}}, Ioss::IJK_t{{57, 9, 65}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A30", zones.back()->m_zone, "zone12", 12, Ioss::IJK_t{{1, -2, -3}}, Ioss::IJK_t{{1, 73, 33}},
      Ioss::IJK_t{{57, 81, 33}}, Ioss::IJK_t{{1, 9, 65}}, Ioss::IJK_t{{57, 1, 65}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "8x48x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B4", zones.back()->m_zone, "zone01", 1, Ioss::IJK_t{{3, 2, -1}}, Ioss::IJK_t{{1, 49, 33}},
      Ioss::IJK_t{{9, 49, 97}}, Ioss::IJK_t{{129, 1, 65}}, Ioss::IJK_t{{65, 1, 73}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B7", zones.back()->m_zone, "zone02", 2, Ioss::IJK_t{{2, -3, -1}}, Ioss::IJK_t{{1, 49, 1}},
      Ioss::IJK_t{{9, 49, 33}}, Ioss::IJK_t{{33, 65, 65}}, Ioss::IJK_t{{1, 73, 65}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B10", zones.back()->m_zone, "zone03", 3, Ioss::IJK_t{{3, -2, 1}}, Ioss::IJK_t{{1, 49, 97}},
      Ioss::IJK_t{{9, 49, 129}}, Ioss::IJK_t{{1, 65, 65}}, Ioss::IJK_t{{33, 65, 73}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B20", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{-1, 2, -3}}, Ioss::IJK_t{{9, 1, 1}},
      Ioss::IJK_t{{9, 49, 129}}, Ioss::IJK_t{{57, 129, 129}}, Ioss::IJK_t{{57, 177, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B25", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{1, 1, 97}},
      Ioss::IJK_t{{9, 1, 129}}, Ioss::IJK_t{{57, 1, 33}}, Ioss::IJK_t{{57, 9, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B27", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{-2, 1, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{9, 1, 33}}, Ioss::IJK_t{{57, 81, 1}}, Ioss::IJK_t{{57, 73, 33}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A31", zones.back()->m_zone, "zone10", 10, Ioss::IJK_t{{-2, -3, 1}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 17, 129}}, Ioss::IJK_t{{1, 1, 17}}, Ioss::IJK_t{{129, 1, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A32", zones.back()->m_zone, "zone10", 10, Ioss::IJK_t{{-3, 2, 1}}, Ioss::IJK_t{{1, 17, 1}},
      Ioss::IJK_t{{1, 33, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{129, 17, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A33", zones.back()->m_zone, "zone10", 10, Ioss::IJK_t{{2, 3, 1}}, Ioss::IJK_t{{1, 33, 1}},
      Ioss::IJK_t{{1, 49, 129}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{129, 17, 17}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "A34", zones.back()->m_zone, "zone12", 12, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{1, 1, 33}},
      Ioss::IJK_t{{9, 1, 97}}, Ioss::IJK_t{{57, 1, 65}}, Ioss::IJK_t{{57, 9, 1}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "128x16x16"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B22", zones.back()->m_zone, "zone07", 7, Ioss::IJK_t{{-2, 1, 3}}, Ioss::IJK_t{{1, 1, 17}},
      Ioss::IJK_t{{129, 17, 17}}, Ioss::IJK_t{{57, 129, 1}}, Ioss::IJK_t{{73, 1, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B31", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{3, -1, -2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{129, 1, 17}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{1, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B32", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{3, 2, -1}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{129, 17, 1}}, Ioss::IJK_t{{1, 17, 1}}, Ioss::IJK_t{{1, 33, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B33", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{3, 1, 2}}, Ioss::IJK_t{{1, 17, 1}},
      Ioss::IJK_t{{129, 17, 17}}, Ioss::IJK_t{{1, 33, 1}}, Ioss::IJK_t{{1, 49, 129}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "32x64x128"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B16", zones.back()->m_zone, "zone05", 5, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 129}},
      Ioss::IJK_t{{33, 65, 129}}, Ioss::IJK_t{{1, 1, 1}}, Ioss::IJK_t{{33, 65, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B17", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{3, -1, -2}}, Ioss::IJK_t{{1, 65, 1}},
      Ioss::IJK_t{{33, 65, 129}}, Ioss::IJK_t{{57, 129, 97}}, Ioss::IJK_t{{57, 1, 129}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B18", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{1, 3, -2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{1, 65, 129}}, Ioss::IJK_t{{57, 129, 33}}, Ioss::IJK_t{{57, 1, 97}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B19", zones.back()->m_zone, "zone06", 6, Ioss::IJK_t{{-3, 1, -2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{33, 1, 129}}, Ioss::IJK_t{{57, 129, 33}}, Ioss::IJK_t{{57, 1, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B26", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{-3, 2, 1}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{33, 65, 1}}, Ioss::IJK_t{{57, 9, 33}}, Ioss::IJK_t{{57, 73, 1}});

  zones.push_back(new Iocgns::StructuredZoneData(zone++, "56x8x64"));
  zones.back()->m_zoneConnectivity.emplace_back(
      "B23", zones.back()->m_zone, "zone07", 7, Ioss::IJK_t{{1, -3, 2}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{57, 1, 65}}, Ioss::IJK_t{{1, 33, 1}}, Ioss::IJK_t{{57, 97, 1}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B28", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{1, 2, 3}}, Ioss::IJK_t{{1, 1, 1}},
      Ioss::IJK_t{{57, 9, 1}}, Ioss::IJK_t{{1, 1, 33}}, Ioss::IJK_t{{57, 9, 33}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B29", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{1, -3, 2}}, Ioss::IJK_t{{1, 9, 1}},
      Ioss::IJK_t{{57, 9, 65}}, Ioss::IJK_t{{1, 9, 33}}, Ioss::IJK_t{{57, 73, 33}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B30", zones.back()->m_zone, "zone08", 8, Ioss::IJK_t{{1, -2, -3}}, Ioss::IJK_t{{1, 1, 65}},
      Ioss::IJK_t{{57, 9, 65}}, Ioss::IJK_t{{1, 81, 33}}, Ioss::IJK_t{{57, 73, 33}});
  zones.back()->m_zoneConnectivity.emplace_back(
      "B34", zones.back()->m_zone, "zone09", 9, Ioss::IJK_t{{2, 1, -3}}, Ioss::IJK_t{{57, 1, 1}},
      Ioss::IJK_t{{57, 9, 65}}, Ioss::IJK_t{{1, 1, 97}}, Ioss::IJK_t{{9, 1, 33}});

  double load_balance_tolerance = 1.2;

  for (size_t proc_count = 3; proc_count <= 384; proc_count *= 2) {
    std::string name = "H1_MK21_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count, 0.75, 1.1); }
  }
  cleanup(zones);
}

TEST_CASE("bc-257x129x2")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  // Failing for line decomposition on 84 processors; 72 works
  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "257x129x2"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;

  double load_balance_tolerance = 1.2;

  for (size_t proc_count = 4; proc_count <= 84; proc_count += 4) {
    std::string name = "BC_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count, 0.9, 1.1); }
  }
  cleanup(zones);
}

TEST_CASE("carnes-mesh")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  // Failing for decomposition on 64 processors
  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "66x2x200"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::J;

  double load_balance_tolerance = 1.2;

  for (size_t proc_count = 2; proc_count <= 64; proc_count *= 2) {
    std::string name = "Carnes_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count); }
  }
  cleanup(zones);
}

TEST_CASE("carnes-blunt-wedge")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  // Failing for decomposition on 64 processors
  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "80x74x1"));
  zones.back()->m_lineOrdinal = Iocgns::Ordinal::K;

  double load_balance_tolerance = 1.2;

  for (size_t proc_count = 2; proc_count <= 64; proc_count *= 2) {
    std::string name = "Carnes_BW_ProcCount_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count); }
  }
  cleanup(zones);
}

TEST_CASE("64GiElem")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  int zone = 1;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "4096x4096x4096"));

  double load_balance_tolerance = 1.01;

  for (size_t proc_count = 2; proc_count <= 1 << 15; proc_count *= 2) {
    std::string name = "64GiElem_PC_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count); }
  }
  cleanup(zones);
}

TEST_CASE("LotsOfZones")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  for (int zone = 1; zone <= 10000; zone++) {
    zones.push_back(new Iocgns::StructuredZoneData(zone, "1000x1000x100"));
  }

  double load_balance_tolerance = 1.01;

  for (size_t proc_count = 2; proc_count <= 1024; proc_count *= 4) {
    std::string name = "Lots_PC_" + std::to_string(proc_count);
    SECTION(name) { check_split_assign(zones, load_balance_tolerance, proc_count, 0.9, 1.1); }
  }
  cleanup(zones);
}

TEST_CASE("half_sphere")
{
  int                                       zone = 1;
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "80x50x24"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "80x50x24"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "80x50x24"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "80x50x50"));
  zones.push_back(new Iocgns::StructuredZoneData(zone++, "24x80x50"));

  double load_balance_tolerance = 1.4;

  std::string name = "half_sphere_8";
  SECTION(name) { check_split_assign(zones, load_balance_tolerance, 8, 0.9, 1.1); }
  cleanup(zones);
}
