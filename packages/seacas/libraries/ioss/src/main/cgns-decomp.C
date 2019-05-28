#include <Ionit_Initializer.h>
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#undef NDEBUG
#include <Ioss_CodeTypes.h>
#include <Ioss_DatabaseIO.h>
#include <Ioss_GetLongOpt.h>
#include <Ioss_IOFactory.h>
#include <Ioss_Property.h>
#include <Ioss_Region.h>
#include <Ioss_ScopeGuard.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_Utils.h>
#include <Ioss_ZoneConnectivity.h>

#include <cgns/Iocgns_StructuredZoneData.h>
#include <cgns/Iocgns_Utils.h>

#include <fmt/format.h>

namespace {
  class Interface
  {
  public:
    bool parse_options(int argc, char **argv)
    {
      int option_index = options_.parse(argc, argv);
      if (option_index < 1) {
        return false;
      }

      if (options_.retrieve("help") != nullptr) {
        options_.usage(std::cerr);
        exit(EXIT_SUCCESS);
      }

      {
        const char *temp = options_.retrieve("processors");
        if (temp != nullptr) {
          proc_count = std::stoi(temp);
        }
        else {
          fmt::print(
              "\nERROR: Processor count must be specified with '-processors $<val>' option.\n");
          options_.usage(std::cerr);
          exit(EXIT_FAILURE);
        }
      }

      {
        const char *temp = options_.retrieve("ordinal");
        if (temp != nullptr) {
          ordinal = std::stoi(temp);
          if (ordinal < 1 || ordinal > 2) {
            fmt::print("\nERROR: Invalid ordinal specified ({}). Must be 0, 1, or 2.\n", ordinal);
            exit(EXIT_FAILURE);
          }
        }
      }

      {
        const char *temp = options_.retrieve("load_balance");
        if (temp != nullptr) {
          load_balance = std::stod(temp);
        }
      }

      if (option_index < argc) {
        filename = argv[option_index];
      }
      else {
        fmt::print(stderr, "\nERROR: Filename missing.\n");
        options_.usage(std::cerr);
        return false;
      }

      return true;
    }

    Interface()
    {
      options_.usage("[options] input_file");
      options_.enroll("help", Ioss::GetLongOption::NoValue, "Print this summary and exit", nullptr);
      options_.enroll("processors", Ioss::GetLongOption::MandatoryValue, "Number of processors.",
                      nullptr);
      options_.enroll("ordinal", Ioss::GetLongOption::MandatoryValue,
                      "Ordinal not to split (0,1,2).", nullptr);
      options_.enroll("load_balance", Ioss::GetLongOption::MandatoryValue,
                      "Max ratio of processor work to average.", nullptr);
    }
    Ioss::GetLongOption options_;
    int                 proc_count{0};
    int                 ordinal{-1};
    double              load_balance{1.4};
    std::string         filename;
  };
} // namespace
namespace {
  std::string codename;
  std::string version = "0.9";

  void cleanup(std::vector<Iocgns::StructuredZoneData *> &zones)
  {
    for (auto &zone : zones) {
      delete zone;
      zone = nullptr;
    }
  }

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

  void decompose(std::vector<Iocgns::StructuredZoneData *> &zones, double load_balance_tolerance,
                 size_t proc_count)

  {
    size_t proc_width = std::floor(std::log10(proc_count)) + 1;
    double total_work =
        std::accumulate(zones.begin(), zones.end(), 0.0,
                        [](double a, Iocgns::StructuredZoneData *b) { return a + b->work(); });

    size_t work_width = std::floor(std::log10(total_work)) + 1;
    work_width += (work_width - 1) / 3; // for the commas...

    double avg_work = total_work / (double)proc_count;
    {
      size_t zcount   = zones.size();
      double max_work = avg_work * load_balance_tolerance;

      Iocgns::Utils::pre_split(zones, avg_work, load_balance_tolerance, 0, proc_count);

      for (const auto zone : zones) {
        if (zone->is_active()) {
          if (zone->work() > max_work) {
            fmt::print("WARNING: Work on zone {} ({:n}) exceeds max_work ({:n}).\n", zone->m_name,
                       zone->work(), (size_t)max_work);
          }
        }
      }

      for (size_t i = 0; i < zones.size(); i++) {
        SMART_ASSERT(zones[i]->m_zone == (int)i + 1)(zones[i]->m_zone)(i + 1);
      }

      {
        std::vector<size_t> work_vector(proc_count);
        Iocgns::Utils::assign_zones_to_procs(zones, work_vector);

        // Print work/processor map...
        std::vector<size_t> proc_work(proc_count);

        fmt::print("\nDecomposition for {} zones over {} processors; Total work = {:n}, Average = "
                   "{:n}\n",
                   zcount, proc_count, (size_t)total_work, (size_t)avg_work);

        // Get max name length for all zones...
        size_t name_len = 0;
        for (const auto zone : zones) {
          if (zone->is_active()) {
            auto len = zone->m_name.length();
            if (len > name_len) {
              name_len = len;
            }
          }
        }

        //=======================================================================
        fmt::print("\n");
        for (const auto adam_zone : zones) {
          if (adam_zone->m_parent == nullptr) {
            if (adam_zone->m_child1 == nullptr) {
              // Unsplit...
              fmt::print("\tZone {:{}}\t  Proc: {:{}}\tOrdinal: {:^12}\tWork: {:{}n} (unsplit)\n",
                         adam_zone->m_name, name_len, adam_zone->m_proc, proc_width,
                         fmt::format("{}x{}x{}", adam_zone->m_ordinal[0], adam_zone->m_ordinal[1],
                                     adam_zone->m_ordinal[2]),
                         adam_zone->work(), work_width);
            }
            else {
              fmt::print("\tZone: {:{}} is decomposed. \tOrdinal: {:^12}\tWork: {:{}n}\n",
                         adam_zone->m_name, name_len,
                         fmt::format("{}x{}x{}", adam_zone->m_ordinal[0], adam_zone->m_ordinal[1],
                                     adam_zone->m_ordinal[2]),
                         adam_zone->work(), work_width);
              for (const auto zone : zones) {
                if (zone->is_active() && zone->m_adam == adam_zone) {
                  fmt::print("\t      {:{}}\t  Proc: {:{}}\tOrdinal: {:^12}\tWork: {:{}n}\n",
                             zone->m_name, name_len, zone->m_proc, proc_width,
                             fmt::format("{}x{}x{}", zone->m_ordinal[0], zone->m_ordinal[1],
                                         zone->m_ordinal[2]),
                             zone->work(), work_width);
                }
              }
            }
          }
        }
        //=======================================================================

        for (const auto zone : zones) {
          if (zone->is_active()) {
            proc_work[zone->m_proc] += zone->work();
          }
        }

        auto v1 = *std::min_element(proc_work.begin(), proc_work.end());
        auto v2 = *std::max_element(proc_work.begin(), proc_work.end());
        if (v1 == v2) {
          fmt::print("\nWork on all processors is {:n}\n\n", v1);
        }
        else {
          int max_star = 40;
          int min_star = max_star * ((double)v1 / (double)(v2));
          int delta    = max_star - min_star;

          fmt::print("\nWork per processor: Minimum = {:n}, Maximum = {:n}, Ratio = {}\n\n", v1, v2,
                     (double)(v2) / v1);
          for (size_t i = 0; i < proc_work.size(); i++) {
            int         star_cnt = (double)(proc_work[i] - v1) / (v2 - v1) * delta + min_star;
            std::string stars(star_cnt, '*');
            fmt::print("\tProcessor {:{}}, work = {:{}n}\t{}\n", i, proc_width, proc_work[i],
                       work_width, stars);
          }
        }

        // Validate decomposition...

        // Each active zone must be on a processor
        for (const auto zone : zones) {
          if (zone->is_active()) {
            SMART_ASSERT(zone->m_proc >= 0)(zone->m_proc);
          }
        }

        // A processor cannot have more than one zone with the same adam zone
        std::set<std::pair<int, int>> proc_adam_map;
        for (const auto zone : zones) {
          if (zone->is_active()) {
            auto success = proc_adam_map.insert(std::make_pair(zone->m_adam->m_zone, zone->m_proc));
            SMART_ASSERT(success.second);
          }
        }

        // Zone Grid Connectivity Checks:
        update_zgc_data(zones, proc_count);

        // Zone Grid Connectivity instances can't connect to themselves...
        for (auto &zone : zones) {
          if (zone->is_active()) {
            for (const auto &zgc : zone->m_zoneConnectivity) {
              if (zgc.is_active()) {
                SMART_ASSERT(zgc.m_ownerZone != zgc.m_donorZone)(zgc.m_ownerZone)(zgc.m_donorZone);
                SMART_ASSERT(zgc.m_ownerGUID != zgc.m_donorGUID)(zgc.m_ownerGUID)(zgc.m_donorGUID);
              }
            }
          }
        }

        // In Iocgns::Utils::common_write_meta_data, there is code to make
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
              SMART_ASSERT(kk.second < 26 * 26 + 26 + 1)(kk.second);
            }
          }
        } //

        // Zone Grid Connectivity from_decomp instances must be symmetric...
        // The GUID encodes the id and the processor,
        std::map<std::pair<size_t, size_t>, int> is_symm;
        for (auto &zone : zones) {
          if (zone->is_active()) {
            for (const auto &zgc : zone->m_zoneConnectivity) {
              if (zgc.is_active() && zgc.is_from_decomp()) {
                is_symm[std::make_pair(std::min(zgc.m_ownerGUID, zgc.m_donorGUID),
                                       std::max(zgc.m_ownerGUID, zgc.m_donorGUID))]++;
              }
            }
          }
        }
        // Iterate `is_symm` and make sure all entries == 2
        for (const auto &item : is_symm) {
          SMART_ASSERT(item.second == 2);
        }
      }
    }
  }
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

  Interface interface;
  bool      success = interface.parse_options(argc, argv);
  if (!success) {
    exit(EXIT_FAILURE);
  }

  std::string in_type = "cgns";

  codename   = argv[0];
  size_t ind = codename.find_last_of('/', codename.size());
  if (ind != std::string::npos) {
    codename = codename.substr(ind + 1, codename.size());
  }

  Ioss::Init::Initializer io;

  Ioss::PropertyManager properties;
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(in_type, interface.filename, Ioss::READ_RESTART,
                                                  (MPI_Comm)MPI_COMM_WORLD, properties);
  if (dbi == nullptr || !dbi->ok()) {
    fmt::print("\nERROR: Could not open database '{}' of type '{}'\n", interface.filename, in_type);
    std::exit(EXIT_FAILURE);
  }

  // NOTE: 'region' owns 'db' pointer at this time...
  Ioss::Region region(dbi, "region_1");

  // Get the structured blocks...
  const auto &blocks = region.get_structured_blocks();
  if (blocks.empty()) {
    fmt::print("\nERROR: There are no structured blocks on the mesh.\n");
    return EXIT_FAILURE;
  }

  std::vector<Iocgns::StructuredZoneData *> zones;
  for (auto iblock : blocks) {
    size_t ni   = iblock->get_property("ni").get_int();
    size_t nj   = iblock->get_property("nj").get_int();
    size_t nk   = iblock->get_property("nk").get_int();
    size_t zone = iblock->get_property("zone").get_int();

    zones.push_back(new Iocgns::StructuredZoneData(iblock->name(), zone, ni, nj, nk));
    if (interface.ordinal >= 0) {
      zones.back()->m_lineOrdinal = interface.ordinal;
    }
  }

  region.output_summary(std::cout, false);
  decompose(zones, interface.load_balance, interface.proc_count);
  cleanup(zones);
}
