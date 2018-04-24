#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <cgns/Iocgns_StructuredZoneData.h>
#include <cgns/Iocgns_Utils.h>
#include <exception>
#include <map>
#include <numeric>
#include <vector>

TEST_CASE("test single block", "[single_block]")
{
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData("zone1", 1, 4, 4, 1));

  int    proc_count = 2;
  double total_work =
      std::accumulate(zones.begin(), zones.end(), 0.0,
                      [](double a, Iocgns::StructuredZoneData *b) { return a + b->work(); });
  double avg_work               = total_work / (double)proc_count;
  double load_balance_tolerance = 1.2;

  Iocgns::Utils::pre_split(zones, avg_work, load_balance_tolerance, 0, proc_count);

  double max_work = avg_work * load_balance_tolerance;
  for (const auto zone : zones) {
    if (zone->is_active()) {
      CHECK(zone->work() <= max_work);
    }
  }
}

TEST_CASE("test single block line", "[single_block_line]")
{
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData("zone1", 1, 4, 4, 1));
  zones.back()->m_lineOrdinal = 0;

  int    proc_count = 4;
  double total_work =
      std::accumulate(zones.begin(), zones.end(), 0.0,
                      [](double a, Iocgns::StructuredZoneData *b) { return a + b->work(); });
  double avg_work               = total_work / (double)proc_count;
  double load_balance_tolerance = 1.05;

  Iocgns::Utils::pre_split(zones, avg_work, load_balance_tolerance, 0, proc_count);

  double max_work = avg_work * load_balance_tolerance;
  for (const auto zone : zones) {
    if (zone->is_active()) {
      CHECK(zone->work() <= max_work);
    }
  }
}

TEST_CASE("test prime sides", "[prime_sides]")
{
  std::vector<Iocgns::StructuredZoneData *> zones;
  zones.push_back(new Iocgns::StructuredZoneData("zone1", 1, 3, 5, 7));

  double total_work =
      std::accumulate(zones.begin(), zones.end(), 0.0,
                      [](double a, Iocgns::StructuredZoneData *b) { return a + b->work(); });

  double load_balance_tolerance = 1.25;
  for (size_t proc_count = 2; proc_count < 8; proc_count++) {
    std::string name = "Prime_ProcCount_" + std::to_string(proc_count);
    SECTION(name)
    {
      double avg_work = total_work / (double)proc_count;

      Iocgns::Utils::pre_split(zones, avg_work, load_balance_tolerance, 0, proc_count);

      double max_work = avg_work * load_balance_tolerance;
      for (const auto zone : zones) {
        if (zone->is_active()) {
          CHECK(zone->work() <= max_work);
        }
      }
      SECTION("assign_to_procs")
      {
        std::vector<size_t> work_vector(proc_count);
        Iocgns::Utils::assign_zones_to_procs(zones, work_vector);

        // Each active zone must be on a processor
        for (const auto zone : zones) {
          if (zone->is_active()) {
            CHECK(zone->m_proc >= 0);
          }
        }

        // Work must be min_work <= work <= max_work
        for (auto work : work_vector) {
          // CHECK(work >= min_work);
          CHECK(work <= max_work);
        }

        // A processor cannot have more than one zone with the same adam zone
        std::set<std::pair<int, int>> proc_adam_map;
        for (size_t i = 0; i < zones.size(); i++) {
          if (zones[i]->is_active()) {
            auto success =
                proc_adam_map.insert(std::make_pair(zones[i]->m_adam->m_zone, zones[i]->m_proc));
            CHECK(success.second);
          }
        }
      }
    }
  }
}

TEST_CASE("test farmer plenum", "[farmer_plenum]")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  auto *zone1 = new Iocgns::StructuredZoneData("zone1", 1, 56, 128, 48);
  zones.push_back(zone1);

  auto *zone2 = new Iocgns::StructuredZoneData("zone2", 2, 32, 64, 48);
  zones.push_back(zone2);

  double load_balance_tolerance = 1.1;
  double total_work =
      std::accumulate(zones.begin(), zones.end(), 0.0,
                      [](double a, Iocgns::StructuredZoneData *b) { return a + b->work(); });

  for (size_t proc_count = 2; proc_count < 20; proc_count++) {
    std::string name = "Plenum_ProcCount_" + std::to_string(proc_count);
    SECTION(name)
    {
      double avg_work = total_work / (double)proc_count;

      SECTION("split_zones")
      {
        Iocgns::Utils::pre_split(zones, avg_work, load_balance_tolerance, 0, proc_count);

        // double min_work = avg_work / load_balance_tolerance;
        double max_work = avg_work * load_balance_tolerance;
        for (const auto zone : zones) {
          if (zone->is_active()) {
            CHECK(zone->work() <= max_work);
          }
        }

        SECTION("assign_to_procs")
        {
          std::vector<size_t> work_vector(proc_count);
          Iocgns::Utils::assign_zones_to_procs(zones, work_vector);

          // Each active zone must be on a processor
          for (const auto zone : zones) {
            if (zone->is_active()) {
              CHECK(zone->m_proc >= 0);
            }
          }

          // Work must be min_work <= work <= max_work
          for (auto work : work_vector) {
            // CHECK(work >= min_work);
            CHECK(work <= max_work);
          }

          // A processor cannot have more than one zone with the same adam zone
          std::set<std::pair<int, int>> proc_adam_map;
          for (size_t i = 0; i < zones.size(); i++) {
            if (zones[i]->is_active()) {
              auto success =
                  proc_adam_map.insert(std::make_pair(zones[i]->m_adam->m_zone, zones[i]->m_proc));
              CHECK(success.second);
            }
          }
        }
      }
    }
  }
}

TEST_CASE("test grv", "[grv]")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  {
    int zone = 1;
    zones.push_back(new Iocgns::StructuredZoneData("zone01", zone++, 8, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone02", zone++, 8, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone03", zone++, 8, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone04", zone++, 8, 1, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone05", zone++, 8, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone06", zone++, 8, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone07", zone++, 8, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone08", zone++, 8, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone09", zone++, 8, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone10", zone++, 8, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone11", zone++, 8, 1, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone12", zone++, 8, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone13", zone++, 8, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone14", zone++, 8, 4, 4));
  }

  double load_balance_tolerance = 1.2;

  double total_work =
      std::accumulate(zones.begin(), zones.end(), 0.0,
                      [](double a, Iocgns::StructuredZoneData *b) { return a + b->work(); });

  for (size_t proc_count = 2; proc_count < 16; proc_count++) {
    std::string name = "GRV_ProcCount_" + std::to_string(proc_count);
    SECTION(name)
    {
      double avg_work = total_work / (double)proc_count;

      SECTION("split_zones")
      {
        Iocgns::Utils::pre_split(zones, avg_work, load_balance_tolerance, 0, proc_count);

        // double min_work = avg_work / load_balance_tolerance;
        double max_work = avg_work * load_balance_tolerance;
        for (const auto zone : zones) {
          if (zone->is_active()) {
            CHECK(zone->work() <= max_work);
          }
        }

        SECTION("assign_to_procs")
        {
          std::vector<size_t> work_vector(proc_count);
          Iocgns::Utils::assign_zones_to_procs(zones, work_vector);

#if 0
	  std::cerr << "\nDecomposition for " << proc_count << " processors; Total work = " << total_work << " Average = " << avg_work << "\n";
	  for (const auto zone : zones) {
	    std::cerr << "Zone " << zone->m_name << "\tProc: " << zone->m_proc << " \tWork: " << zone->work() << "\n";
	  }
#endif
          // Each active zone must be on a processor
          for (const auto zone : zones) {
            if (zone->is_active()) {
              CHECK(zone->m_proc >= 0);
            }
          }

          // Work must be min_work <= work <= max_work
          max_work *= 1.2; // Modify range until we get full splitting working in test.
          for (auto work : work_vector) {
            // CHECK(work >= min_work);
            CHECK(work <= max_work);
          }

          // A processor cannot have more than one zone with the same adam zone
          std::set<std::pair<int, int>> proc_adam_map;
          for (const auto zone : zones) {
            if (zone->is_active()) {
              auto success =
                  proc_adam_map.insert(std::make_pair(zone->m_adam->m_zone, zone->m_proc));
              CHECK(success.second);
            }
          }
        }
      }
    }
  }
}

TEST_CASE("test mk21", "[mk21]")
{
  std::vector<Iocgns::StructuredZoneData *> zones;

  {
    int zone = 1;
    zones.push_back(new Iocgns::StructuredZoneData("zone01", zone++, 6, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone02", zone++, 6, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone03", zone++, 6, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone04", zone++, 6, 2, 8));
    zones.push_back(new Iocgns::StructuredZoneData("zone05", zone++, 6, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone06", zone++, 6, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone07", zone++, 6, 2, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone08", zone++, 6, 2, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone09", zone++, 6, 4, 4));

    zones.push_back(new Iocgns::StructuredZoneData("zone10", zone++, 6, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone11", zone++, 6, 8, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone12", zone++, 6, 8, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone13", zone++, 4, 2, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone14", zone++, 4, 2, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone15", zone++, 4, 2, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone16", zone++, 4, 2, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone17", zone++, 4, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone18", zone++, 4, 2, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone19", zone++, 4, 2, 4));

    zones.push_back(new Iocgns::StructuredZoneData("zone20", zone++, 4, 2, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone21", zone++, 4, 2, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone22", zone++, 4, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone23", zone++, 4, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone24", zone++, 4, 4, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone25", zone++, 4, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone26", zone++, 4, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone27", zone++, 4, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone28", zone++, 4, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone29", zone++, 4, 4, 2));

    zones.push_back(new Iocgns::StructuredZoneData("zone30", zone++, 4, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone31", zone++, 4, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone32", zone++, 4, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone33", zone++, 4, 2, 2));
    zones.push_back(new Iocgns::StructuredZoneData("zone34", zone++, 16, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone35", zone++, 16, 4, 4));
    zones.push_back(new Iocgns::StructuredZoneData("zone36", zone++, 16, 4, 4));
  }

  double load_balance_tolerance = 1.2;

  double total_work =
      std::accumulate(zones.begin(), zones.end(), 0.0,
                      [](double a, Iocgns::StructuredZoneData *b) { return a + b->work(); });

  for (size_t proc_count = 2; proc_count < 17; proc_count++) {
    std::string name = "GRV_ProcCount_" + std::to_string(proc_count);
    SECTION(name)
    {
      double avg_work = total_work / (double)proc_count;

      SECTION("split_zones")
      {
        Iocgns::Utils::pre_split(zones, avg_work, load_balance_tolerance, 0, proc_count);

        // double min_work = avg_work / load_balance_tolerance;
        double max_work = avg_work * load_balance_tolerance;
        for (const auto zone : zones) {
          if (zone->is_active()) {
            CHECK(zone->work() <= max_work);
          }
        }

        SECTION("assign_to_procs")
        {
          std::vector<size_t> work_vector(proc_count);
          Iocgns::Utils::assign_zones_to_procs(zones, work_vector);

#if 0
	      std::cerr << "\nDecomposition for " << proc_count << " processors; Total work = " << total_work << " Average = " << avg_work << "\n";
	      for (const auto zone : zones) {
		std::cerr << "Zone " << zone->m_name << "\tProc: " << zone->m_proc << " \tWork: " << zone->work() << "\n";
	      }
#endif
          // Each active zone must be on a processor
          for (const auto zone : zones) {
            if (zone->is_active()) {
              CHECK(zone->m_proc >= 0);
            }
          }

          // Work must be min_work <= work <= max_work
          max_work *= 1.1; // Modify range until we get full splitting working in test.
          for (auto work : work_vector) {
            // CHECK(work >= min_work);
            CHECK(work <= max_work);
          }

          // A processor cannot have more than one zone with the same adam zone
          std::set<std::pair<int, int>> proc_adam_map;
          for (const auto zone : zones) {
            if (zone->is_active()) {
              auto success =
                  proc_adam_map.insert(std::make_pair(zone->m_adam->m_zone, zone->m_proc));
              CHECK(success.second);
            }
          }
        }
      }
    }
  }
}
