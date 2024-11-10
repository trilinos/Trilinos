// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#include <cstdlib>
#include <exception>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "add_to_log.h"
#define FMT_DEPRECATED_OSTREAM
#include "fmt/ostream.h"
#include "fmt/ranges.h"
#include "format_time.h"
#include "hwm.h"
#include "open_file_limit.h"
#include "time_stamp.h"

#include <Ionit_Initializer.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>

#include <cgns/Iocgns_Utils.h>

#include "CP_SystemInterface.h"
#include "CP_Version.h"

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif

unsigned int debug_level = 0;

#if FMT_VERSION >= 90000
#include "Ioss_ZoneConnectivity.h"
#include "Ioss_StructuredBlock.h"

namespace fmt {
  template <> struct formatter<Ioss::ZoneConnectivity> : ostream_formatter
  {
  };
} // namespace fmt
namespace fmt {
  template <> struct formatter<Ioss::BoundaryCondition> : ostream_formatter
  {
  };
} // namespace fmt
#endif

namespace {
  std::string tsFormat = "[{:%H:%M:%S}] ";

  using GlobalZgcMap   = std::map<std::pair<std::string, std::string>, Ioss::ZoneConnectivity>;
  using GlobalBcMap    = std::map<std::pair<std::string, std::string>, Ioss::BoundaryCondition>;
  using GlobalBlockMap = std::map<const std::string, const Ioss::StructuredBlock *>;
  using GlobalIJKMap   = std::map<const std::string, Ioss::IJK_t>;
  using PartVector     = std::vector<std::unique_ptr<Ioss::Region>>;

  GlobalZgcMap generate_global_zgc(const PartVector &part_mesh);
  GlobalBcMap  generate_global_bc(const PartVector &part_mesh);

  void   info_structuredblock(const Ioss::Region &region);
  void   resolve_offsets(const PartVector &part_mesh, GlobalBlockMap &all_blocks);
  void   update_global_ijk(const PartVector &part_mesh, GlobalIJKMap &global_block);
  void   transfer_nodal_field(const Ioss::StructuredBlock *sb, const std::vector<double> &input,
                              std::vector<double> &output);
  void   transfer_cell_field(const Ioss::StructuredBlock *sb, const std::vector<double> &input,
                             std::vector<double> &output);
  void   transfer_nodal_coordinates(const PartVector &part_mesh, Ioss::Region &output_region,
                                    bool minimize_open_files);
  double transfer_step(const PartVector &part_mesh, Ioss::Region &output_region, int istep,
                       bool minimize_open_files);
  void   union_zgc_range(Ioss::ZoneConnectivity &zgc_i, const Ioss::ZoneConnectivity &zgc_j);
  void   union_bc_range(Ioss::IJK_t &g_beg, Ioss::IJK_t &g_end, const Ioss::IJK_t &l_beg,
                        const Ioss::IJK_t &l_end, const Ioss::IJK_t &offset);

  int get_constant_face(const Ioss::IJK_t &beg, const Ioss::IJK_t &end)
  {
    for (int i = 0; i < 3; i++) {
      if (beg[i] == end[i]) {
        return (beg[i] == 1) ? i : i + 3;
      }
    }
    return 7;
  }

  bool is_field_valid(const Cpup::StringVector &variable_list, const std::string &field_name)
  {
    if (variable_list.empty() ||
        (variable_list.size() == 1 && Ioss::Utils::str_equal(variable_list[0], "all"))) {
      return true;
    }

    // At this point, the variable_list contains one or more entries
    // of fields that should be output on combined file.  Run through
    // list and see if `field_name` is in the list.
    return std::any_of(
        variable_list.begin(), variable_list.end(),
        [&field_name](const auto &valid) { return Ioss::Utils::str_equal(valid, field_name); });
  }

  int verify_timestep_count(const PartVector &part_mesh)
  {
    int num_time_steps = part_mesh[0]->get_property("state_count").get_int();

    bool differ = false;
    for (const auto &part : part_mesh) {
      int nts = part->get_property("state_count").get_int();
      if (nts != num_time_steps) {
        differ = true;
      }
      num_time_steps = num_time_steps < nts ? num_time_steps : nts;
    }
    if (differ) {
      fmt::print(stderr,
                 "\nWARNING: The number of time steps is not the same on all input databases.\n"
                 "         Using minimum count of {}\n\n",
                 num_time_steps);
    }
    else {
      fmt::print(stderr, "\nNumber of time steps on input databases = {}\n\n", num_time_steps);
    }
    return num_time_steps;
  }

} // namespace

template <typename INT> void cpup(Cpup::SystemInterface &interFace, INT dummy);

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  try {
    Cpup::SystemInterface::show_version();
    Ioss::Init::Initializer io;

    Cpup::SystemInterface interFace;
    bool                  ok = interFace.parse_options(argc, argv);

    debug_level = interFace.debug();

    if (!ok) {
      fmt::print(stderr, "\nERROR: Problems parsing command line arguments.\n\n");
      exit(EXIT_FAILURE);
    }

    int    error = 0;
    double begin = Ioss::Utils::timer();
    cpup(interFace, static_cast<int64_t>(0));
    double end = Ioss::Utils::timer();

    fmt::print(stderr,
               "\nTotal Execution Time = {:.2f} seconds, Maximum memory = {} MiBytes.\n******* "
               "END *******\n\n",
               end - begin,
               fmt::group_digits((get_hwm_memory_info() + 1024 * 1024 - 1) / (1024 * 1024)));

    add_to_log(argv[0], end - begin);

#ifdef SEACAS_HAVE_MPI
    MPI_Finalize();
#endif

    return (error);
  }
  catch (std::exception &e) {
    fmt::print(stderr, "ERROR: Standard exception: {}\n", e.what());
  }
  exit(EXIT_SUCCESS);
}

template <typename INT> void cpup(Cpup::SystemInterface &interFace, INT /*dummy*/)
{
  auto width = Ioss::Utils::number_width(interFace.processor_count(), false);

  bool minimize_open_files = interFace.minimize_open_files();

  if (!minimize_open_files) {
    // Query the system to see if the number of files exceeds the system limit and we
    // need to force use of minimize_open_files...
    int max_files = open_file_limit() - 1; // We also have an output file.
    if (interFace.processor_count() > max_files) {
      minimize_open_files = true;
      fmt::print("Single file mode... (Max open = {})\n", max_files);
    }
  }

  PartVector part_mesh(interFace.processor_count());
  for (int p = 0; p < interFace.processor_count(); p++) {
    std::string root_dir = interFace.root_dir();
    std::string sub_dir  = interFace.sub_dir();
    std::string prepend{};

    if (!root_dir.empty()) {
      prepend = root_dir + "/";
    }
    else if (Ioss::Utils::is_path_absolute(prepend)) {
      prepend = "";
    }
    else {
      prepend = "./";
    }
    if (!sub_dir.empty()) {
      prepend += sub_dir + "/";
    }

    prepend += interFace.basename() + "." + interFace.cgns_suffix();
    auto filename = Ioss::Utils::decode_filename(prepend, p, interFace.processor_count());

    if (debug_level & 1) {
      fmt::print(stderr, "{} Processor rank {:{}}, file {}\n", time_stamp(tsFormat), p, width,
                 filename);
    }
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create("cgns", filename, Ioss::READ_RESTART,
                                                    Ioss::ParallelUtils::comm_world());
    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    dbi->set_field_separator(1);

    // NOTE: region owns database pointer at this time...
    std::string name = "CPUP_" + std::to_string(p + 1);
    part_mesh[p]     = std::make_unique<Ioss::Region>(dbi, name);

    if (part_mesh[p]->mesh_type() != Ioss::MeshType::STRUCTURED) {
      part_mesh[p]->output_summary(std::cerr);
      fmt::print(stderr,
                 "\nERROR: Only UNSTRUCTURED CGNS file joining is supported at this time.\n");
      exit(EXIT_FAILURE);
    }

    if (debug_level & 2) {
      part_mesh[p]->output_summary(std::cerr);
      fmt::print(stderr, "\n");
    }
    if (minimize_open_files) {
      part_mesh[p]->get_database()->closeDatabase();
    }
  }

  // Each processor may have a different set of zones.  This routine
  // will sync the information such that at end, there is a consistent
  // set of structuredBlocks defined with the correct local and
  // global, i,j,k ranges and offsets.

  GlobalBlockMap all_blocks;
  GlobalIJKMap   global_block;
  for (const auto &part : part_mesh) {
    const auto &blocks = part->get_structured_blocks();
    for (const auto &block : blocks) {
      const auto &name = block->name();
      all_blocks[name] = block;

      // Build map of unique blocks in the mesh.
      auto [part_name, proc] = Iocgns::Utils::decompose_name(name, true);
      global_block[part_name];

      // Set the zgc 'from_decomp' property...
      for (const auto &zgc : block->m_zoneConnectivity) {
        const auto [zgc_name, zgc_proc] = Iocgns::Utils::decompose_name(zgc.m_donorName, true);
        if (zgc_name == part_name) {
          zgc.m_fromDecomp = true;
        }
      }
    }
  }

  // Resolve the offsets for each split block based on the decomp-ZGC
  // with the same block on other processors
  resolve_offsets(part_mesh, all_blocks);

  // Offsets are now correct. Need to calculate and update the global size of each unique block...
  update_global_ijk(part_mesh, global_block);

  // Need a consistent set of Boundary Conditions for each zone unioned across the proc-local-zones
  GlobalBcMap global_bc = generate_global_bc(part_mesh);

  // Need a consistent set of ZGC for each zone unioned across the proc-local-zones...
  // Skip the ZGC that are "from_decomp"
  GlobalZgcMap global_zgc = generate_global_zgc(part_mesh);

  // Create output file...
  Ioss::PropertyManager properties{};
  properties.add(Ioss::Property("FLUSH_INTERVAL", 0));
  Ioss::DatabaseIO *dbo =
      Ioss::IOFactory::create("cgns", interFace.output_filename(), Ioss::WRITE_RESTART,
                              Ioss::ParallelUtils::comm_world(), properties);
  if (dbo == nullptr || !dbo->ok(true)) {
    std::exit(EXIT_FAILURE);
  }

  // NOTE: 'output_region' owns 'dbo' pointer at this time
  Ioss::Region output_region(dbo, "cpup_output_region");
  output_region.property_add(Ioss::Property("code_name", qainfo[0]));
  output_region.property_add(Ioss::Property("code_version", qainfo[1] + ":" + qainfo[2]));

  output_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

  // KLUGE: Remove this...
  // Doesn't affect the output model at all, so not a big issue, ...
  auto *nb = new Ioss::NodeBlock(dbo, "nodeblock_1", 1, 3);
  output_region.add(nb);

  // Create the output structured blocks...
  for (auto &block_range : global_block) {
    const auto &block_name = block_range.first;
    auto       *block      = new Ioss::StructuredBlock(dbo, block_name, 3, block_range.second);
    output_region.add(block);

    // Add BC to the block...
    for (const auto &bc_map : global_bc) {
      if (bc_map.first.first == block_name) {
        block->m_boundaryConditions.push_back(bc_map.second);
      }
    }

    // Add ZGC to the block...
    for (const auto &zgc_map : global_zgc) {
      if (zgc_map.first.first == block_name) {
        block->m_zoneConnectivity.push_back(zgc_map.second);
      }
    }
  }

  // Copy the sidesets and assemblies from the proc-0 input file to the output file...
  auto       &part  = part_mesh[0];
  const auto &ssets = part->get_sidesets();
  for (const auto &sset : ssets) {
    auto *oss = new Ioss::SideSet(*sset);
    output_region.add(oss);
  }

  const auto &assems = part->get_assemblies();
  for (const auto &assem : assems) {
    auto *oass = new Ioss::Assembly(*assem);
    output_region.add(oass);
  }

  if (debug_level & 4) {
    info_structuredblock(output_region);
  }

  output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

  output_region.begin_mode(Ioss::STATE_MODEL);
  transfer_nodal_coordinates(part_mesh, output_region, minimize_open_files);
  output_region.end_mode(Ioss::STATE_MODEL);

  // ******* Transient Data...
  output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

  // ... Iterate the output_region structured blocks,
  //     .. Find corresponding structured blocks on part meshes.
  //        .. Add each valid block and node_block field
  const auto &variable_list = interFace.var_names();
  if (!(variable_list.size() == 1 && Ioss::Utils::str_equal(variable_list[0], "none") == 0)) {
    const auto &blocks = output_region.get_structured_blocks();
    for (const auto &block : blocks) {
      int64_t num_cell = block->get_property("cell_count").get_int();
      int64_t num_node = block->get_property("node_count").get_int();

      auto &onb = block->get_node_block();

      // Find all corresponding blocks on the input part meshes...
      for (const auto &prt : part_mesh) {
        const auto &pblocks = prt->get_structured_blocks();
        for (const auto &pblock : pblocks) {
          const auto &name            = pblock->name();
          auto [part_name, part_proc] = Iocgns::Utils::decompose_name(name, true);
          if (part_name == block->name()) {
            Ioss::NameList fields = pblock->field_describe(Ioss::Field::TRANSIENT);

            for (const auto &field_name : fields) {
              if (is_field_valid(variable_list, field_name)) {
                Ioss::Field field = pblock->get_field(field_name);
                if (!block->field_exists(field_name)) {
                  // If the field does not already exist, add it to the output block...
                  field.reset_count(num_cell);
                  block->field_add(field);
                }
              }
            }

            // Now transfer the fields on the embedded node block...
            auto pnb = pblock->get_node_block();
            fields.clear();
            pnb.field_describe(Ioss::Field::TRANSIENT, &fields);
            for (const auto &field_name : fields) {
              if (is_field_valid(variable_list, field_name)) {
                Ioss::Field field = pnb.get_field(field_name);
                if (!onb.field_exists(field_name)) {
                  field.reset_count(num_node);
                  // If the field does not already exist, add it to the output block...
                  onb.field_add(field);
                }
              }
            }

            break; // Should be only a single instance of each block on a part mesh.
          }
        }
      }
    }
  }

  output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

  output_region.begin_mode(Ioss::STATE_TRANSIENT);

  int num_time_steps = verify_timestep_count(part_mesh);

  // Determine if user wants a subset of timesteps transferred to the output file.
  int ts_min  = interFace.step_min();
  int ts_max  = interFace.step_max();
  int ts_step = interFace.step_interval();

  if (ts_min < 0) {
    ts_min = num_time_steps + 1 + ts_min;
  }
  if (ts_max < 0) {
    ts_max = num_time_steps + 1 + ts_max;
  }

  // Time steps for output file
  int time_step_out = 0;

  ts_max = ts_max < num_time_steps ? ts_max : num_time_steps;
  if (ts_min <= ts_max) {
    fmt::print(stderr, "\tTransferring step {} to step {} by {}\n", ts_min, ts_max, ts_step);
  }

  // Determine how many steps will be written...
  int output_steps = (ts_max - ts_min) / ts_step + 1;

  double start_time = Ioss::Utils::timer();
  double cur_time   = start_time;
  for (int time_step = ts_min; time_step <= ts_max; time_step += ts_step) {
    time_step_out++;
    double time_val = transfer_step(part_mesh, output_region, time_step, minimize_open_files);

    double time_per_step       = Ioss::Utils::timer() - cur_time;
    cur_time                   = Ioss::Utils::timer();
    double elapsed             = cur_time - start_time;
    double avg_time_per_step   = elapsed / time_step_out;
    double percentage_done     = (time_step_out * 100.0) / output_steps;
    double estimated_remaining = avg_time_per_step * (output_steps - time_step_out);
    if (debug_level & 1) {
      fmt::print(stderr,
                 "{} \tWrote step {:6}, time {:8.4e}\t[{:5.1f}%, Elapsed={}, \tETA={}, \tTPS={}]\n",
                 time_stamp(tsFormat), fmt::group_digits(time_step), time_val, percentage_done,
                 format_time(elapsed), format_time(estimated_remaining),
                 format_time(time_per_step));
    }
    else {
      fmt::print(stderr,
                 "\tWrote step {:6}, time {:8.4e}\t[{:5.1f}%, Elapsed={}, ETA={}, TPS={}]       \r",
                 fmt::group_digits(time_step), time_val, percentage_done, format_time(elapsed),
                 format_time(estimated_remaining), format_time(time_per_step));
    }
  }
  output_region.end_mode(Ioss::STATE_TRANSIENT);

  fmt::print(stderr, "\n\n********************* OUTPUT DATABASE ********************\n");
  output_region.output_summary(std::cerr);
}

namespace {
  GlobalZgcMap generate_global_zgc(const PartVector &part_mesh)
  {
    GlobalZgcMap global_zgc;
    for (const auto &part : part_mesh) {
      const auto &blocks = part->get_structured_blocks();
      for (const auto &block : blocks) {
        const auto [part_name, proc] = Iocgns::Utils::decompose_name(block->name(), true);

        for (const auto &zgc : block->m_zoneConnectivity) {
          if (!zgc.m_fromDecomp) {
            auto &gzgc = global_zgc[std::make_pair(part_name, zgc.m_connectionName)];
            if (gzgc.m_connectionName.empty()) {
              // First time this ZGC has been found.  Copy from the per-proc instance and update...
              gzgc.m_connectionName = zgc.m_connectionName;
              gzgc.m_donorName      = Iocgns::Utils::decompose_name(zgc.m_donorName, true).first;
              gzgc.m_transform      = zgc.m_transform;
            }

            // Create a temporary zgc; adjust its ranges to "global" by adding the offsets
            // and then union if with the global instance...
            Ioss::IJK_t own_off = block->get_ijk_offset();

            auto tmp_zgc{zgc};
            tmp_zgc.m_ownerRangeBeg[0] += own_off[0];
            tmp_zgc.m_ownerRangeBeg[1] += own_off[1];
            tmp_zgc.m_ownerRangeBeg[2] += own_off[2];
            tmp_zgc.m_ownerRangeEnd[0] += own_off[0];
            tmp_zgc.m_ownerRangeEnd[1] += own_off[1];
            tmp_zgc.m_ownerRangeEnd[2] += own_off[2];

            // Now find the donor block...
            auto  donor       = Iocgns::Utils::decompose_name(zgc.m_donorName, true);
            auto *donor_block = part_mesh[donor.second]->get_structured_block(zgc.m_donorName);
            SMART_ASSERT(donor_block != nullptr);
            Ioss::IJK_t don_off = donor_block->get_ijk_offset();

            tmp_zgc.m_donorRangeBeg[0] += don_off[0];
            tmp_zgc.m_donorRangeBeg[1] += don_off[1];
            tmp_zgc.m_donorRangeBeg[2] += don_off[2];
            tmp_zgc.m_donorRangeEnd[0] += don_off[0];
            tmp_zgc.m_donorRangeEnd[1] += don_off[1];
            tmp_zgc.m_donorRangeEnd[2] += don_off[2];

            union_zgc_range(gzgc, tmp_zgc);
          }
        }
      }
    }
    return global_zgc;
  }

  GlobalBcMap generate_global_bc(const PartVector &part_mesh)
  {
    GlobalBcMap global_bc;
    for (const auto &part : part_mesh) {
      const auto &blocks = part->get_structured_blocks();
      for (const auto &block : blocks) {
        Ioss::IJK_t offset           = block->get_ijk_offset();
        const auto [part_name, proc] = Iocgns::Utils::decompose_name(block->name(), true);
        const auto &sb_bc            = block->m_boundaryConditions;
        for (const auto &bc : sb_bc) {
          auto &gbc = global_bc[std::make_pair(part_name, bc.m_bcName)];
          if (gbc.m_bcName.empty()) {
            gbc.m_bcName  = bc.m_bcName;
            gbc.m_famName = bc.m_famName;
          }
          union_bc_range(gbc.m_rangeBeg, gbc.m_rangeEnd, bc.m_rangeBeg, bc.m_rangeEnd, offset);
          SMART_ASSERT(gbc.which_face() == bc.which_face());
        }
      }
    }
    return global_bc;
  }

  double transfer_step(const PartVector &part_mesh, Ioss::Region &output_region, int istep,
                       bool minimize_open_files)
  {
    double time  = part_mesh[0]->get_state_time(istep);
    int    ostep = output_region.add_state(time);

    output_region.begin_state(ostep);

    for (const auto &part : part_mesh) {
      part->begin_state(istep);
    }

    const auto &blocks = output_region.get_structured_blocks();
    for (const auto &block : blocks) {
      int64_t             num_cell = block->get_property("cell_count").get_int();
      std::vector<double> output(num_cell);
      std::vector<double> input;

      Ioss::NameList fields = block->field_describe(Ioss::Field::TRANSIENT);

      // Not sure if this is the best ordering of loops, but it minimizes the
      // amount of data gathered at one time at the cost of multiple iterations
      // through the block-finding loop...
      for (const auto &field_name : fields) {

        // Find all corresponding blocks on the input part meshes...
        for (const auto &part : part_mesh) {
          const auto &pblocks = part->get_structured_blocks();
          for (const auto &pblock : pblocks) {
            const auto &name             = pblock->name();
            const auto [part_name, proc] = Iocgns::Utils::decompose_name(name, true);
            if (part_name == block->name()) {
              if (pblock->field_exists(field_name)) {
                pblock->get_field_data(field_name, input);
                transfer_cell_field(pblock, input, output);
              }
              break; // Should be only a single instance of each block on a part mesh.
            }
          }
          if (minimize_open_files) {
            part->get_database()->closeDatabase();
          }
        }
        block->put_field_data(field_name, output);
      }
    }

    // Now do the fields on the embedded node block...
    for (const auto &block : blocks) {
      int64_t num_node = block->get_property("node_count").get_int();
      auto   &onb      = block->get_node_block();

      std::vector<double> output(num_node);
      std::vector<double> input;

      Ioss::NameList fields = onb.field_describe(Ioss::Field::TRANSIENT);

      // Not sure if this is the best ordering of loops, but it minimizes the
      // amount of data gathered at one time at the cost of multiple iterations
      // through the block-finding loop...
      for (const auto &field_name : fields) {

        // Find all corresponding blocks on the input part meshes...
        for (const auto &part : part_mesh) {
          const auto &pblocks = part->get_structured_blocks();
          for (const auto &pblock : pblocks) {
            const auto &name             = pblock->name();
            const auto [part_name, proc] = Iocgns::Utils::decompose_name(name, true);
            if (part_name == block->name()) {
              auto &inb = pblock->get_node_block();
              if (inb.field_exists(field_name)) {
                inb.get_field_data(field_name, input);
                transfer_nodal_field(pblock, input, output);
              }
              break; // Should be only a single instance of each block on a part mesh.
            }
          }
          if (minimize_open_files) {
            part->get_database()->closeDatabase();
          }
        }
        onb.put_field_data(field_name, output);
      }
    }

    output_region.end_state(ostep);
    for (const auto &part : part_mesh) {
      part->end_state(istep);
    }
    return time;
  }

  void resolve_offsets(const PartVector &part_mesh, GlobalBlockMap &all_blocks)
  {
    bool change_made;
    do {
      change_made = false;
      for (const auto &part : part_mesh) {
        const auto &blocks = part->get_structured_blocks();
        for (const auto &block : blocks) {
          for (const auto &zgc : block->m_zoneConnectivity) {
            if (zgc.is_from_decomp()) {
              auto plane = get_constant_face(zgc.m_ownerRangeBeg, zgc.m_ownerRangeEnd);
              if (plane < 3) {
                // This zone connects to another zone "below" it.
                // Find the connecting zone and adjust the correct offset
                const auto *donor     = all_blocks[zgc.m_donorName];
                auto        offset    = donor->get_ijk_offset();
                auto        range     = donor->get_ijk_local();
                auto        my_offset = block->get_ijk_offset();
                if (my_offset[plane] != range[plane] + offset[plane]) {
                  block->set_ijk_offset(plane, range[plane] + offset[plane]);
                  change_made = true;
                }
              }
            }
          }
        }
      }
    } while (change_made);
  }

  void update_global_ijk(const PartVector &part_mesh, GlobalIJKMap &global_block)
  {
    for (const auto &part : part_mesh) {
      const auto &blocks = part->get_structured_blocks();
      for (const auto &block : blocks) {
        auto ijk_o                   = block->get_ijk_offset();
        auto ijk_g                   = block->get_ijk_global();
        const auto [part_name, proc] = Iocgns::Utils::decompose_name(block->name(), true);
        auto &cur_global             = global_block[part_name];
        cur_global[0]                = std::max(ijk_o[0] + ijk_g[0], cur_global[0]);
        cur_global[1]                = std::max(ijk_o[1] + ijk_g[1], cur_global[1]);
        cur_global[2]                = std::max(ijk_o[2] + ijk_g[2], cur_global[2]);
      }
    }

    for (const auto &part : part_mesh) {
      const auto &blocks = part->get_structured_blocks();
      for (const auto &block : blocks) {
        const auto [part_name, proc] = Iocgns::Utils::decompose_name(block->name(), true);
        auto &cur_global             = global_block[part_name];
        block->set_ijk_global(cur_global);
      }
      if (debug_level & 4) {
        info_structuredblock(*part);
      }
    }
  }

  void info_structuredblock(const Ioss::Region &region)
  {
    const Ioss::StructuredBlockContainer &sbs = region.get_structured_blocks();
    for (const auto &sb : sbs) {

      auto ijk_global = sb->get_ijk_global();
      fmt::print(stderr, "\n{} '{}' {}", sb->type_string(), sb->name(), fmt::join(ijk_global, "x"));

      auto ijk_offset = sb->get_ijk_offset();
      auto ijk_local  = sb->get_ijk_local();
      fmt::print(stderr, " [{}, Offset = {}] ", fmt::join(ijk_local, "x"),
                 fmt::join(ijk_offset, ", "));

      int64_t num_cell = sb->get_property("cell_count").get_int();
      int64_t num_node = sb->get_property("node_count").get_int();
      fmt::print(stderr, "  {:14} cells, {:14} nodes ", fmt::group_digits(num_cell),
                 fmt::group_digits(num_node));

#if defined(__NVCC__)
#define CONST
#else
#define CONST const
#endif
      if (!sb->m_zoneConnectivity.empty()) {
        fmt::print(stderr, "\n\tConnectivity with other blocks:\n");
        for (CONST auto &zgc : sb->m_zoneConnectivity) {
          fmt::print(stderr, "{}\n", zgc);
        }
      }

      if (!sb->m_boundaryConditions.empty()) {
        fmt::print(stderr, "\tBoundary Conditions:\n");
        for (CONST auto &bc : sb->m_boundaryConditions) {
          fmt::print(stderr, "{}\n", bc);
        }
      }
#undef CONST
    }
  }

  void transfer_nodal_field(const Ioss::StructuredBlock *sb, const std::vector<double> &input,
                            std::vector<double> &output)
  {
    // Get the IJK subset that this nodeblock covers...
    auto ijkg = sb->get_ijk_global();
    auto ijkl = sb->get_ijk_local();

    // Get the IJK offset...
    auto offset = sb->get_ijk_offset();

    SMART_ASSERT(input.size() == (size_t)(ijkl[0] + 1) * (ijkl[1] + 1) * (ijkl[2] + 1))
    (input.size())(ijkl[0])(ijkl[1])(ijkl[2]);
    SMART_ASSERT(output.size() == (size_t)(ijkg[0] + 1) * (ijkg[1] + 1) * (ijkg[2] + 1))
    (output.size())(ijkg[0])(ijkg[1])(ijkg[2]);

    size_t idx = 0;
    for (int k = offset[2]; k < offset[2] + ijkl[2] + 1; k++) {
      for (int j = offset[1]; j < offset[1] + ijkl[1] + 1; j++) {
        for (int i = offset[0]; i < offset[0] + ijkl[0] + 1; i++) {
          size_t oidx  = i + (ijkg[0] + 1) * j + (ijkg[0] + 1) * (ijkg[1] + 1) * k;
          output[oidx] = input[idx++];
        }
      }
    }
  }

  void transfer_cell_field(const Ioss::StructuredBlock *sb, const std::vector<double> &input,
                           std::vector<double> &output)
  {
    // Get the IJK subset that this structured covers...
    auto ijkg = sb->get_ijk_global();
    auto ijkl = sb->get_ijk_local();

    // Get the IJK offset...
    auto offset = sb->get_ijk_offset();

    SMART_ASSERT(input.size() == (size_t)ijkl[0] * ijkl[1] * ijkl[2])
    (input.size())(ijkl[0])(ijkl[1])(ijkl[2]);
    SMART_ASSERT(output.size() == (size_t)ijkg[0] * ijkg[1] * ijkg[2])
    (output.size())(ijkg[0])(ijkg[1])(ijkg[2]);

    size_t idx = 0;
    for (int k = offset[2]; k < offset[2] + ijkl[2]; k++) {
      for (int j = offset[1]; j < offset[1] + ijkl[1]; j++) {
        for (int i = offset[0]; i < offset[0] + ijkl[0]; i++) {
          size_t oidx  = i + ijkg[0] * j + ijkg[0] * ijkg[1] * k;
          output[oidx] = input[idx++];
        }
      }
    }
  }

  void transfer_nodal_coordinates(const PartVector &part_mesh, Ioss::Region &output_region,
                                  bool minimize_open_files)
  {
    // This implementation results in having to iterate over the part
    // mesh 3 times -- once for each coordinate axis, but minimizes
    // the memory requirements.  Will monitor whether this results in
    // excessive computation time, but does permit larger models to be
    // processed.

    std::array<std::string, 3> fields{"mesh_model_coordinates_x", "mesh_model_coordinates_y",
                                      "mesh_model_coordinates_z"};

    const auto &blocks = output_region.get_structured_blocks();
    for (const auto &block : blocks) {
      // Get size of node_block...
      auto               &onb       = block->get_node_block();
      size_t              num_coord = onb.entity_count();
      std::vector<double> coord(num_coord);
      for (int dim = 0; dim < 3; dim++) {

        // Find all corresponding blocks on the input part meshes...
        for (const auto &part : part_mesh) {
          const auto &pblocks = part->get_structured_blocks();
          for (const auto &pblock : pblocks) {
            const auto &name             = pblock->name();
            const auto [part_name, proc] = Iocgns::Utils::decompose_name(name, true);
            if (part_name == block->name()) {
              std::vector<double> lcoord;
              pblock->get_field_data(fields[dim], lcoord);
              transfer_nodal_field(pblock, lcoord, coord);
              break; // Should be only a single instance of each block on a part mesh.
            }
          }
          if (minimize_open_files) {
            part->get_database()->closeDatabase();
          }
        }
        block->put_field_data(fields[dim], coord);
      }
    }
  }

  void union_bc_range(Ioss::IJK_t &g_beg, Ioss::IJK_t &g_end, const Ioss::IJK_t &l_beg,
                      const Ioss::IJK_t &l_end, const Ioss::IJK_t &offset)
  {
    for (int i = 0; i < 3; i++) {
      g_beg[i] = g_beg[i] == 0 ? l_beg[i] + offset[i] : std::min(g_beg[i], l_beg[i] + offset[i]);
      g_end[i] = std::max(g_end[i], l_end[i] + offset[i]);
    }
  }

  void union_zgc_range(Ioss::ZoneConnectivity &zgc_i, const Ioss::ZoneConnectivity &zgc_j)
  {
    assert(zgc_i.m_transform == zgc_j.m_transform);
    if (zgc_i.m_ownerRangeBeg[0] == 0 && zgc_i.m_ownerRangeBeg[1] == 0 &&
        zgc_i.m_ownerRangeBeg[2] == 0) {
      // This is a newly-created zgc that hasn't been unioned with anything yet.
      zgc_i.m_ownerRangeBeg = zgc_j.m_ownerRangeBeg;
      zgc_i.m_ownerRangeEnd = zgc_j.m_ownerRangeEnd;
      zgc_i.m_donorRangeBeg = zgc_j.m_donorRangeBeg;
      zgc_i.m_donorRangeEnd = zgc_j.m_donorRangeEnd;
    }
    else {
      for (int i = 0; i < 3; i++) {
        if (zgc_i.m_ownerRangeBeg[i] <= zgc_i.m_ownerRangeEnd[i]) {
          zgc_i.m_ownerRangeBeg[i] = std::min(zgc_i.m_ownerRangeBeg[i], zgc_j.m_ownerRangeBeg[i]);
          zgc_i.m_ownerRangeEnd[i] = std::max(zgc_i.m_ownerRangeEnd[i], zgc_j.m_ownerRangeEnd[i]);
        }
        else {
          zgc_i.m_ownerRangeBeg[i] = std::max(zgc_i.m_ownerRangeBeg[i], zgc_j.m_ownerRangeBeg[i]);
          zgc_i.m_ownerRangeEnd[i] = std::min(zgc_i.m_ownerRangeEnd[i], zgc_j.m_ownerRangeEnd[i]);
        }

        if (zgc_i.m_donorRangeBeg[i] <= zgc_i.m_donorRangeEnd[i]) {
          zgc_i.m_donorRangeBeg[i] = std::min(zgc_i.m_donorRangeBeg[i], zgc_j.m_donorRangeBeg[i]);
          zgc_i.m_donorRangeEnd[i] = std::max(zgc_i.m_donorRangeEnd[i], zgc_j.m_donorRangeEnd[i]);
        }
        else {
          zgc_i.m_donorRangeBeg[i] = std::max(zgc_i.m_donorRangeBeg[i], zgc_j.m_donorRangeBeg[i]);
          zgc_i.m_donorRangeEnd[i] = std::min(zgc_i.m_donorRangeEnd[i], zgc_j.m_donorRangeEnd[i]);
        }
      }
    }
  }
} // namespace
