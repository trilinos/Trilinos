// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

// Make asserts active even in non-debug build
#undef NDEBUG

#include <Ionit_Initializer.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
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
#include <sys/ioctl.h>
#include <unistd.h>
#include <utility>
#include <vector>

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

#include <fmt/color.h>
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

      verbose = options_.retrieve("verbose") != nullptr;

      if (options_.retrieve("output") != nullptr) {
        const std::string temp = options_.retrieve("output");
        histogram              = temp.find("h") != std::string::npos;
        work_per_processor     = temp.find("w") != std::string::npos;
        zone_proc_assignment   = temp.find("z") != std::string::npos;
        verbose                = temp.find("v") != std::string::npos || verbose;
        communication_map      = temp.find("c") != std::string::npos;
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
        const char *temp = options_.retrieve("line_decomposition");
        if (temp != nullptr) {
          line_decomposition = temp;
        }
      }

      {
        const char *temp = options_.retrieve("db_type");
        if (temp != nullptr) {
          filetype = temp;
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
        options_.usage(std::cout);
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
      options_.enroll("line_decomposition", Ioss::GetLongOption::MandatoryValue,
                      "list of 1 or more BC (Family) names.\n"
                      "\t\tFor all structured zones which this BC touches, the ordinal of the face "
                      "(i,j,k) will\n"
                      "\t\tbe set such that a parallel decomposition will not split the zone along "
                      "this ordinal.",
                      nullptr);
      options_.enroll("load_balance", Ioss::GetLongOption::MandatoryValue,
                      "Max ratio of processor work to average. [default 1.4]", nullptr);
      options_.enroll("verbose", Ioss::GetLongOption::NoValue,
                      "Print additional decomposition information", nullptr);
      options_.enroll("db_type", Ioss::GetLongOption::MandatoryValue,
                      "Database Type: gen_struc or cgns. Default is cgns.", nullptr);
      options_.enroll("output", Ioss::GetLongOption::MandatoryValue,
                      "What is printed: z=zone-proc assignment, h=histogram, w=work-per-processor, "
                      "c=comm map, v=verbose.",
                      "zhwc");
    }
    Ioss::GetLongOption options_;
    int                 proc_count{0};
    int                 ordinal{-1};
    double              load_balance{1.4};
    std::string         filename{};
    std::string         filetype{"cgns"};
    std::string         line_decomposition{};
    bool                verbose{false};
    bool                histogram{true};
    bool                work_per_processor{true};
    bool                zone_proc_assignment{true};
    bool                communication_map{true};
  };
} // namespace
namespace {
  std::string codename;
  std::string version = "0.96";

  int term_width();

  void cleanup(std::vector<Iocgns::StructuredZoneData *> &zones)
  {
    for (auto &zone : zones) {
      delete zone;
      zone = nullptr;
    }
  }

  double surface_ratio(const Iocgns::StructuredZoneData *zone)
  {
    size_t surf =
        2 * (zone->m_ordinal[0] * zone->m_ordinal[1] + zone->m_ordinal[0] * zone->m_ordinal[2] +
             zone->m_ordinal[1] * zone->m_ordinal[2]);
    size_t vol = zone->cell_count();

    // If a 'perfect' cube, then would be pl=cbrt(vol) on a side and surf would be 6*pl*pl
    // Calculate 'perfect' surf / actual surf...
    double pl    = std::cbrt(vol);
    double psurf = 6.0 * pl * pl;
    return double(surf) / psurf;
  }

  int64_t generate_guid(size_t id, int rank, int proc_count)
  {
    static size_t lpow2 = 0;
    static int    procs = -1;
    if (procs != proc_count) {
      lpow2 = Ioss::Utils::log_power_2(proc_count);
      procs = proc_count;
    }
    SMART_ASSERT(rank >= 0);
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

  void validate_decomposition(std::vector<Iocgns::StructuredZoneData *> &zones, int proc_count)
  {

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

  void output_communications(std::vector<Iocgns::StructuredZoneData *> &zones, int proc_count)
  {
    fmt::print("Communication Map: [] is from decomposition; () is from zone-to-zone; omits "
               "on-proc communication.\n");
    for (const auto &adam_zone : zones) {
      if (adam_zone->m_parent == nullptr) {
        std::vector<std::pair<int, int>> comms;

        // Iterate children (or self) of the adam_zone.
        for (const auto zone : zones) {
          if (zone->is_active() && zone->m_adam == adam_zone) {
            for (auto &zgc : zone->m_zoneConnectivity) {
              int p1 = zgc.m_ownerProcessor;
              int p2 = zgc.m_donorProcessor;
              if (p1 != p2) {
                int pmin = std::min(p1, p2);
                int pmax = std::max(p1, p2);
                if (zgc.is_from_decomp()) {
                  comms.emplace_back(pmin, -pmax);
                }
                else {
                  comms.emplace_back(pmin, pmax);
                }
              }
            }
          }
        }

        Ioss::Utils::uniquify(comms);

        int pw = Ioss::Utils::number_width(proc_count, false);
        // Two tabs at beginning ~16 spaces.  Each entry is "[pw->pw]  " which is 6+2pw
        int npl  = (term_width() - 16) / (6 + 2 * pw);
        npl      = npl < 1 ? 1 : npl;
        int line = 0;

        fmt::print("\tZone '{}' ({} inter-processor communications):\n\t\t", adam_zone->m_name,
                   comms.size());
        for (const auto &proc : comms) {
          if (proc.second < 0) {
            // From decomposition
            fmt::print(fg(fmt::color::yellow), "[{:{}}->{:{}}]  ", proc.first, pw, -proc.second,
                       pw);
          }
          else {
            // Zone to Zone
            fmt::print("({:{}}->{:{}})  ", proc.first, pw, proc.second, pw);
          }
          if (++line >= npl) {
            fmt::print("\n\t\t");
            line = 0;
          }
        }
        fmt::print("\n");
        if (line > 0) {
          fmt::print("\n");
        }
      }
    }
  }

  void output_histogram(const std::vector<size_t> &proc_work, size_t avg_work, size_t median)
  {
    fmt::print("Work-per-processor Histogram\n");
    std::array<size_t, 16> histogram{};

    auto wmin = *std::min_element(proc_work.begin(), proc_work.end());
    auto wmax = *std::max_element(proc_work.begin(), proc_work.end());

    size_t hist_size = std::min(size_t(16), (wmax - wmin));
    hist_size        = std::min(hist_size, proc_work.size());

    if (hist_size <= 1) {
      fmt::print("\tWork is the same on all processors; no histogram needed.\n\n");
      return;
    }

    auto delta = double(wmax + 1 - wmin) / hist_size;
    for (const auto &pw : proc_work) {
      auto bin = size_t(double(pw - wmin) / delta);
      SMART_ASSERT(bin < hist_size)(bin)(hist_size);
      histogram[bin]++;
    }

    size_t proc_width = Ioss::Utils::number_width(proc_work.size(), true);
    size_t work_width = Ioss::Utils::number_width(wmax, true);

    fmt::print("\n\t{:^{}} {:^{}}\n", "Work Range", 2 * work_width + 2, "#", proc_width);
    auto hist_max = *std::max_element(histogram.begin(), histogram.end());
    for (size_t i = 0; i < hist_size; i++) {
      int         max_star = 50;
      int         star_cnt = ((double)histogram[i] / hist_max * max_star);
      std::string stars(star_cnt, '*');
      for (int j = 9; j < star_cnt;) {
        stars[j] = '|';
        j += 10;
      }
      if (histogram[i] > 0 && star_cnt == 0) {
        stars = '.';
      }
      size_t      w1 = wmin + size_t(i * delta);
      size_t      w2 = wmin + size_t((i + 1) * delta);
      std::string postfix;
      if (w1 <= avg_work && avg_work < w2) {
        postfix += "average";
      }
      if (w1 <= median && median < w2) {
        if (!postfix.empty()) {
          postfix += ", ";
        }
        postfix += "median";
      }
      fmt::print("\t{:{}n}..{:{}n} ({:{}n}):\t{:{}}  {}\n", w1, work_width, w2, work_width,
                 histogram[i], proc_width, stars, max_star, postfix);
    }
    fmt::print("\n");
  }
  void describe_decomposition(std::vector<Iocgns::StructuredZoneData *> &zones,
                              size_t orig_zone_count, const Interface &interFace)
  {
    size_t proc_count = interFace.proc_count;
    bool   verbose    = interFace.verbose;

    // ========================================================================
    // Output information about decomposition...
    double total_work = std::accumulate(zones.begin(), zones.end(), 0.0,
                                        [](double a, Iocgns::StructuredZoneData *b) {
                                          return a + (b->m_parent == nullptr ? b->work() : 0);
                                        });

    // Get some information just to make output look better.  Not part of actual decomposition.
    size_t proc_width = Ioss::Utils::number_width(proc_count, false);
    size_t work_width = Ioss::Utils::number_width((size_t)total_work, true);

    // Find maximum ordinal to get width... (makes output look better)
    int max_ordinal = 0;
    for (const auto zone : zones) {
      if (zone->is_active()) {
        max_ordinal =
            std::max({max_ordinal, zone->m_ordinal[0], zone->m_ordinal[1], zone->m_ordinal[2]});
      }
    }
    size_t ord_width = Ioss::Utils::number_width(max_ordinal, false);
    double avg_work  = total_work / (double)proc_count;

    // Print work/processor map...
    fmt::print("\nDecomposing {:n} zones over {:n} processors; Total work = {:n}; Average = "
               "{:n} (goal)\n",
               orig_zone_count, proc_count, (size_t)total_work, (size_t)avg_work);

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

    if (interFace.zone_proc_assignment) {
      // Output Zone->Processor assignment info
      fmt::print("\n");
      for (const auto adam_zone : zones) {
        if (adam_zone->m_parent == nullptr) {
          if (adam_zone->m_child1 == nullptr) {
            // Unsplit...
            fmt::print("\tZone: {:{}}\t  Proc: {:{}}\tOrd: {:^12}    Work: {:{}n} (unsplit)\n",
                       adam_zone->m_name, name_len, adam_zone->m_proc, proc_width,
                       fmt::format("{1:{0}} x {2:{0}} x {3:{0}}", ord_width,
                                   adam_zone->m_ordinal[0], adam_zone->m_ordinal[1],
                                   adam_zone->m_ordinal[2]),
                       adam_zone->work(), work_width);
          }
          else {
            fmt::print("\tZone: {:{}} is decomposed. \tOrd: {:^12}    Work: {:{}n}\n",
                       adam_zone->m_name, name_len,
                       fmt::format("{1:{0}} x {2:{0}} x {3:{0}}", ord_width,
                                   adam_zone->m_ordinal[0], adam_zone->m_ordinal[1],
                                   adam_zone->m_ordinal[2]),
                       adam_zone->work(), work_width);
            for (const auto zone : zones) {
              if (zone->is_active() && zone->m_adam == adam_zone) {
                fmt::print("\t      {:{}}\t  Proc: {:{}}\tOrd: {:^12}    Work: {:{}n}    SurfExp: "
                           "{:0.3}\n",
                           zone->m_name, name_len, zone->m_proc, proc_width,
                           fmt::format("{1:{0}} x {2:{0}} x {3:{0}}", ord_width, zone->m_ordinal[0],
                                       zone->m_ordinal[1], zone->m_ordinal[2]),
                           zone->work(), work_width, surface_ratio(zone));
              }
            }
          }
        }
      }
    }

    // Output Processor Work information
    std::vector<size_t> proc_work(proc_count);
    for (const auto zone : zones) {
      if (zone->is_active()) {
        proc_work[zone->m_proc] += zone->work();
      }
    }

    auto   min_work = *std::min_element(proc_work.begin(), proc_work.end());
    auto   max_work = *std::max_element(proc_work.begin(), proc_work.end());
    size_t median   = 0;
    {
      auto pw_copy(proc_work);
      std::nth_element(pw_copy.begin(), pw_copy.begin() + pw_copy.size() / 2, pw_copy.end());
      median = pw_copy[pw_copy.size() / 2];
      fmt::print(
          "\nWork per processor:\n\tMinimum = {:n}, Maximum = {:n}, Median = {:n}, Ratio = {:.3}\n\n",
          min_work, max_work, median, (double)(max_work) / min_work);
    }
    if (interFace.work_per_processor) {
      if (min_work == max_work) {
        fmt::print("\nWork on all processors is {:n}\n\n", min_work);
      }
      else {
        int max_star = 40;
        int min_star = max_star * ((double)min_work / (double)(max_work));
        int delta    = max_star - min_star;

        for (size_t i = 0; i < proc_work.size(); i++) {
          int star_cnt =
              (double)(proc_work[i] - min_work) / (max_work - min_work) * delta + min_star;
          std::string stars(star_cnt, '*');
          std::string format = "\tProcessor {:{}}, work = {:{}n}  ({:.2f})\t{}\n";
          if (proc_work[i] == max_work) {
            fmt::print(fg(fmt::color::red), format, i, proc_width, proc_work[i], work_width,
                       proc_work[i] / avg_work, stars);
          }
          else if (proc_work[i] == min_work) {
            fmt::print(fg(fmt::color::green), format, i, proc_width, proc_work[i], work_width,
                       proc_work[i] / avg_work, stars);
          }
          else {
            fmt::print(format, i, proc_width, proc_work[i], work_width, proc_work[i] / avg_work,
                       stars);
          }
          if (verbose) {
            for (const auto zone : zones) {
              if ((size_t)zone->m_proc == i) {
                auto pct = int(100.0 * (double)zone->work() / proc_work[i] + 0.5);
                fmt::print("\t      {:{}} {:{}n}\t{:3}%\t{:^12}\n", zone->m_name, name_len,
                           zone->work(), work_width, pct,
                           fmt::format("{1:{0}} x {2:{0}} x {3:{0}}", ord_width, zone->m_ordinal[0],
                                       zone->m_ordinal[1], zone->m_ordinal[2]));
              }
            }
            fmt::print("\n");
          }
        }
      }
    }

    // Output Histogram...
    if (interFace.histogram) {
      output_histogram(proc_work, (size_t)avg_work, median);
    }

    // Communication Information (proc X communicates with proc Z)
    if (interFace.communication_map) {
      output_communications(zones, proc_count);
    }

    // Calculate "nodal inflation" -- number of new surface nodes created...
    auto nodal_work = std::accumulate(zones.begin(), zones.end(), (size_t)0,
                                      [](size_t a, Iocgns::StructuredZoneData *b) {
                                        return a + (b->m_parent == nullptr ? b->node_count() : 0);
                                      });

    if (nodal_work > 0) {
      auto new_nodal_work = std::accumulate(zones.begin(), zones.end(), (size_t)0,
                                            [](size_t a, Iocgns::StructuredZoneData *b) {
                                              return a + (b->is_active() ? b->node_count() : 0);
                                            });

      auto delta = new_nodal_work - nodal_work;
      fmt::print("Nodal Inflation:\n\tOriginal Node Count = {:n}, Decomposed Node Count = {:n}, "
                 "Created = {:n}, Ratio = {:.2f}\n\n",
                 nodal_work, new_nodal_work, delta, (double)new_nodal_work / nodal_work);
    }

    // Imbalance penalty -- max work / avg work.  If perfect balance, then all processors would have
    // "avg_work" work to do. With current decomposition, every processor has to wait until
    // "max_work" is done.  Penalty = max_work / avg_work.
    fmt::print("Imbalance Penalty:\n\tMaximum Work = {:n}, Average Work = {:n}, Penalty (max/avg) "
               "= {:.2f}\n\n",
               max_work, (size_t)avg_work, (double)max_work / avg_work);
  }
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

  Interface interFace;
  bool      success = interFace.parse_options(argc, argv);
  if (!success) {
    exit(EXIT_FAILURE);
  }

  std::string in_type = interFace.filetype;

  codename   = argv[0];
  size_t ind = codename.find_last_of('/', codename.size());
  if (ind != std::string::npos) {
    codename = codename.substr(ind + 1, codename.size());
  }

  Ioss::Init::Initializer io;

  Ioss::PropertyManager properties{};
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(in_type, interFace.filename, Ioss::READ_RESTART,
                                                  (MPI_Comm)MPI_COMM_WORLD, properties);
  if (dbi == nullptr || !dbi->ok()) {
    fmt::print("\nERROR: Could not open database '{}' of type '{}'\n", interFace.filename, in_type);
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
    if (interFace.ordinal >= 0) {
      zones.back()->m_lineOrdinal = interFace.ordinal;
    }
    zones.back()->m_zoneConnectivity = iblock->m_zoneConnectivity;
  }

  if (in_type == "cgns") {
    Iocgns::Utils::set_line_decomposition(dbi->get_file_pointer(), interFace.line_decomposition,
                                          zones, 0, interFace.verbose);
  }

  region.output_summary(std::cerr, false);

  size_t orig_zone_count = zones.size();
  Iocgns::Utils::decompose_model(zones, interFace.proc_count, 0, interFace.load_balance,
                                 interFace.verbose);
  update_zgc_data(zones, interFace.proc_count);

  describe_decomposition(zones, orig_zone_count, interFace);

  validate_decomposition(zones, interFace.proc_count);

  cleanup(zones);
}

#if defined(_MSC_VER)
#include <io.h>
#define isatty _isatty
#endif

namespace {
  int term_width()
  {
    int cols = 100;
    if (isatty(fileno(stdout))) {
#ifdef TIOCGSIZE
      struct ttysize ts
      {
      };
      ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
      cols = ts.ts_cols;
#elif defined(TIOCGWINSZ)
      struct winsize ts;
      ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
      cols = ts.ws_col;
#endif /* TIOCGSIZE */
    }
    return cols;
  }
} // namespace
