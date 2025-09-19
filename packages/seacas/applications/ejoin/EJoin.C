// Copyright(C) 1999-2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unistd.h>
#include <vector>

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include "add_to_log.h"
#include "time_stamp.h"

#include <exodusII.h>

#include <Ionit_Initializer.h>
#include <Ioss_Enumerate.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Transform.h>
#include <Ioss_Utils.h>
#include <tokenize.h>

#include "EJ_CodeTypes.h"
#include "EJ_SystemInterface.h"
#include "EJ_Version.h"
#include "EJ_mapping.h"
#include "EJ_match.h"
#include "EJ_vector3d.h"

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif

// Globals...
std::map<Ioss::NodeSet *, std::vector<size_t>> nodeset_in_out_map;
IO_map                                         output_input_map;

namespace {
  const std::string tsFormat    = "[%H:%M:%S] ";
  unsigned int      debug_level = 0;

  bool valid_variable(const std::string &variable, size_t id, const StringIdVector &variable_list);
  bool check_variable_mismatch(const std::string &type, const StringIdVector &variable_list,
                               const Ioss::NameList &fields);

  bool define_global_fields(Ioss::Region &output_region, const RegionVector &part_mesh,
                            const StringIdVector &variable_list);
  bool define_nodal_fields(const Ioss::Region &output_region, const RegionVector &part_mesh,
                           const StringIdVector &variable_list, const SystemInterface &interFace);
  bool define_element_fields(const Ioss::Region   &output_region,
                             const StringIdVector &variable_list);
  bool define_nodeset_fields(const Ioss::Region   &output_region,
                             const StringIdVector &variable_list);
  bool define_sideset_fields(const Ioss::Region   &output_region,
                             const StringIdVector &variable_list);
  void define_nodal_nodeset_fields(const Ioss::Region &output_region, const RegionVector &part_mesh,
                                   const StringIdVector  &variable_list,
                                   const SystemInterface &interFace);

  template <typename INT>
  void output_nodeblock(Ioss::Region &output_region, const RegionVector &part_mesh,
                        const std::vector<INT> &local_node_map, std::vector<INT> &global_node_map);
  template <typename INT>
  void output_elementblock(Ioss::Region &output_region, const RegionVector &part_mesh,
                           const std::vector<INT> &local_node_map,
                           const std::vector<INT> &local_element_map, bool ignore_element_ids);
  template <typename INT>
  void output_nodeset(Ioss::Region &output_region, const std::vector<INT> &local_node_map,
                      bool nodes_consolidated);
  template <typename INT>
  void output_sideset(Ioss::Region &output_region, const std::vector<INT> &local_element_map);
  template <typename INT>
  void output_nodal_nodeset(Ioss::Region &output_region, const RegionVector &part_mesh,
                            const SystemInterface  &interFace,
                            const std::vector<INT> &local_node_map);
  template <typename INT>
  void output_transient_state(Ioss::Region &output_region, const RegionVector &part_mesh,
                              double time, const std::vector<INT> &local_node_map,
                              SystemInterface &interFace, bool merged);
  void process_nodeset_omissions(const RegionVector &part_mesh, const Omissions &omit);
  void process_sideset_omissions(const RegionVector &part_mesh, const Omissions &omit);
  void process_assembly_omissions(const RegionVector &part_mesh, const Omissions &omit);

  int count_omissions(Ioss::Region *region)
  {
    int         omitted = 0;
    const auto &blocks  = region->get_element_blocks();
    for (const auto &block : blocks) {
      if (block->property_exists(std::string("omitted"))) {
        omitted++;
      }
    }
    return omitted;
  }

  template <typename T> bool approx_equal(T v1, T v2, T offset)
  {
#if 1
    static const T tolerance = 100.0f * std::numeric_limits<float>::epsilon();
    double         d1        = std::fabs(v1 - v2);
    double         d2        = std::fabs((v1 - offset) + (v2 - offset)) * tolerance;
    return d1 <= d2;
#else
    return (float)v1 == (float)v2;
#endif
  }

  using EntityIdSet = std::set<std::pair<Ioss::EntityType, int64_t>>;
  EntityIdSet id_set;

  void set_id(Ioss::GroupingEntity *old_ge, Ioss::GroupingEntity *new_ge)
  {
    if (old_ge->property_exists("id")) {
      int64_t id      = old_ge->get_property("id").get_int();
      bool    succeed = id_set.insert(std::make_pair(new_ge->type(), id)).second;
      if (succeed) {
        // Add id as the property "id" to ge
        new_ge->property_add(Ioss::Property("id", id));
      }
    }
  }

  void transfer_elementblock(const Ioss::Region &region, Ioss::Region &output_region,
                             bool create_assemblies, bool combine_similar, bool debug);
  void transfer_assembly(const Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_nodesets(const Ioss::Region &region, Ioss::Region &output_region,
                         bool combine_similar, bool debug);
  void transfer_sidesets(const Ioss::Region &region, Ioss::Region &output_region,
                         bool combine_similar, bool debug);
  void create_nodal_nodeset(const Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role, const std::string &prefix = "");

  // Used when combining two or more input set into a single output set (WIP)
  // May reequire other operations to completely combine the entities...
  // This only handles the entity count.
  void reset_entity_count(Ioss::GroupingEntity *ge, int64_t new_entity_count)
  {
    ge->reset_entity_count(new_entity_count);
    auto field_names = ge->field_describe();
    for (const auto &field_name : field_names) {
      const auto &field_ref = ge->get_fieldref(field_name);
      const_cast<Ioss::Field &>(field_ref).reset_count(new_entity_count);
    }
  }

  // Used when combining two or more input set into a single output set (WIP)
  // May reequire other operations to completely combine the entities...
  // This only handles the entity count.
  void add_to_entity_count(Ioss::GroupingEntity *ge, int64_t entity_count_increment)
  {
    int64_t new_count = ge->entity_count() + entity_count_increment;
    reset_entity_count(ge, new_count);
  }

  template <typename INT>
  void check_for_duplicate_nodeset_nodes(Ioss::Region           &output_region,
                                         const std::vector<INT> &local_node_map);
} // namespace

template <typename INT>
double ejoin(SystemInterface &interFace, const RegionVector &part_mesh, INT dummy);

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  try {
    SystemInterface::show_version();
    Ioss::Init::Initializer io;

    SystemInterface interFace;
    bool            ok = interFace.parse_options(argc, argv);

    if (!ok) {
      fmt::print(stderr, "\nERROR: Problems parsing command line arguments.\n\n");
      exit(EXIT_FAILURE);
    }

    debug_level = interFace.debug();

    if ((debug_level & 64) != 0U) {
      ex_opts(EX_VERBOSE | EX_DEBUG);
    }
    else {
      ex_opts(0);
    }

    int error = 0;

    int int_byte_size = 4;
    if (interFace.ints64bit()) {
      int_byte_size = 8;
    }

    const Omissions                &omissions  = interFace.block_omissions();
    const Omissions                &inclusions = interFace.block_inclusions();
    RegionVector                    part_mesh(interFace.inputFiles_.size());
    std::vector<Ioss::DatabaseIO *> dbi(interFace.inputFiles_.size());
    for (size_t p = 0; p < interFace.inputFiles_.size(); p++) {
      dbi[p] = Ioss::IOFactory::create("exodusII", interFace.inputFiles_[p], Ioss::READ_RESTART,
                                       Ioss::ParallelUtils::comm_world());
      if (dbi[p] == nullptr || !dbi[p]->ok(true)) {
        std::exit(EXIT_FAILURE);
      }

      if (dbi[p]->int_byte_size_api() > int_byte_size) {
        int_byte_size = dbi[p]->int_byte_size_api();
      }
      dbi[p]->set_lowercase_database_names(false);
    }

    for (size_t p = 0; p < interFace.inputFiles_.size(); p++) {
      dbi[p]->set_surface_split_type(Ioss::SPLIT_BY_DONT_SPLIT);

      if (int_byte_size == 8) {
        dbi[p]->set_int_byte_size_api(Ioss::USE_INT64_API);
      }

      if (interFace.disable_field_recognition()) {
        dbi[p]->set_field_separator(1);
      }

      if (!omissions[p].empty() || !inclusions[p].empty()) {
        dbi[p]->set_block_omissions(omissions[p], inclusions[p]);
      }

      // Generate a name for the region based on the part number...
      std::string prefix = interFace.block_prefix();
      std::string name   = prefix + std::to_string(p + 1);
      // NOTE: region owns database pointer at this time...
      part_mesh[p] = new Ioss::Region(dbi[p], name);

      int omission_count = count_omissions(part_mesh[p]);
      part_mesh[p]->property_add(Ioss::Property("block_omission_count", omission_count));

      const vector3d &offset    = interFace.offset(p);
      const vector3d &scale     = interFace.scale(p);
      bool            is_offset = offset.x != 0.0 || offset.y != 0.0 || offset.z != 0.0;
      bool            is_scale  = scale.x != 1.0 || scale.y != 1.0 || scale.z != 1.0;
      if (is_offset || is_scale) {
        Ioss::NodeBlock *nb    = part_mesh[p]->get_node_blocks()[0];
        Ioss::Field      coord = nb->get_field("mesh_model_coordinates");
        if (is_scale) {
          auto *transform = Ioss::Transform::create("scale3D");
          SMART_ASSERT(transform != nullptr);
          std::vector<double> values{scale.x, scale.y, scale.z};
          transform->set_properties("scale", values);
          coord.add_transform(transform);
        }
        if (is_offset) {
          auto *transform = Ioss::Transform::create("offset3D");
          SMART_ASSERT(transform != nullptr);
          std::vector<double> values{offset.x, offset.y, offset.z};
          transform->set_properties("offset", values);
          coord.add_transform(transform);
        }
        nb->field_erase("mesh_model_coordinates");
        nb->field_add(coord);
      }
    }

    process_nodeset_omissions(part_mesh, interFace.nodeset_omissions());
    process_sideset_omissions(part_mesh, interFace.sideset_omissions());
    process_assembly_omissions(part_mesh, interFace.assembly_omissions());

    double time = 0.0;

    if (int_byte_size == 4) {
      time = ejoin(interFace, part_mesh, 0);
    }
    else {
      time = ejoin(interFace, part_mesh, static_cast<int64_t>(0));
    }

    for (auto &pm : part_mesh) {
      delete pm;
    }

    add_to_log(argv[0], time);

#ifdef SEACAS_HAVE_MPI
    MPI_Comm parent_comm;
    MPI_Comm_get_parent(&parent_comm);
    if (parent_comm != MPI_COMM_NULL) {
      MPI_Barrier(parent_comm);
    }
    MPI_Finalize();
#endif

    return error;
  }
  catch (std::exception &e) {
    fmt::print(stderr, "ERROR: Standard exception: {}\n", e.what());
  }
}

void process_specified_combines(const RegionVector &part_mesh, const std::string &split,
                                Ioss::EntityType type)
{
  const auto combines = Ioss::tokenize(split, ";");
  for (const auto &combine : combines) {
    size_t end         = combine.find(':');
    auto   output_name = combine.substr(0, end);
    auto   inputs      = combine.substr(end + 1);
    auto   input_names = Ioss::tokenize(inputs, ",");

    for (const auto *part : part_mesh) {
      for (const auto &name : input_names) {
        auto *entity = part->get_entity(name, type);
        if (entity != nullptr) {
          entity->property_add(Ioss::Property(std::string("ejoin_combine_into"), output_name));
        }
      }
    }
  }
}

template <typename INT>
double ejoin(SystemInterface &interFace, const RegionVector &part_mesh, INT /*dummy*/)
{
  double begin      = Ioss::Utils::timer();
  size_t part_count = interFace.inputFiles_.size();
  SMART_ASSERT(part_count == part_mesh.size());

  Ioss::PropertyManager properties{};
  if (sizeof(INT) == 8) {
    properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
    properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
  }

  if (interFace.use_netcdf4()) {
    properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
  }

  if (interFace.compression_level() > 0 || interFace.szip() || interFace.quantize() ||
      interFace.zlib() || interFace.zstd() || interFace.bz2()) {
    properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    properties.add(Ioss::Property("COMPRESSION_LEVEL", interFace.compression_level()));
    properties.add(Ioss::Property("COMPRESSION_SHUFFLE", 1));

    if (interFace.zlib()) {
      properties.add(Ioss::Property("COMPRESSION_METHOD", "zlib"));
    }
    else if (interFace.szip()) {
      properties.add(Ioss::Property("COMPRESSION_METHOD", "szip"));
    }
    else if (interFace.zstd()) {
      properties.add(Ioss::Property("COMPRESSION_METHOD", "zstd"));
    }
    else if (interFace.bz2()) {
      properties.add(Ioss::Property("COMPRESSION_METHOD", "bzip2"));
    }

    if (interFace.quantize()) {
      properties.add(Ioss::Property("COMPRESSION_QUANTIZE_NSD", interFace.quantize_nsd()));
    }
  }

  // Get maximum length of names used on any of the input files...
  int max_name_length = 32;
  for (const auto *mesh : part_mesh) {
    int max_name_used = mesh->get_database()->maximum_symbol_length();
    max_name_length   = std::max(max_name_length, max_name_used);
  }
  properties.add(Ioss::Property("MAXIMUM_NAME_LENGTH", max_name_length));

  properties.add(Ioss::Property("FLUSH_INTERVAL", 0));
  Ioss::DatabaseIO *dbo =
      Ioss::IOFactory::create("exodusII", interFace.outputName_, Ioss::WRITE_RESTART,
                              Ioss::ParallelUtils::comm_world(), properties);
  if (dbo == nullptr || !dbo->ok(true)) {
    std::exit(EXIT_FAILURE);
  }

  // NOTE: 'output_region' owns 'dbo' pointer at this time
  Ioss::Region output_region(dbo, "ejoin_output_region");

  output_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

  output_region.property_add(Ioss::Property("code_name", qainfo[0]));
  output_region.property_add(Ioss::Property("code_version", qainfo[2]));

  if (debug_level & 1) {
    fmt::print(stderr, "{}", time_stamp(tsFormat));
  }

  INT node_offset = 0;
  for (auto &pm : part_mesh) {
    pm->property_add(Ioss::Property("node_offset", node_offset));
    INT local_node_count = pm->get_property("node_count").get_int();
    node_offset += local_node_count;
  }

  INT              node_count = node_offset; // Sum of nodes in part meshes.
  std::vector<INT> local_node_map(node_count);
  std::vector<INT> global_node_map;

  // Need a map from local node to global node.  Position in the map
  // is part_mesh.offset+local_node_position Return value is position
  // in the global node list.  If no mapping, then the list is simply
  // 0..number_of_global_nodes where number_of_global_nodes is the sum
  // of the node counts in the individual files.
  // This routine also eliminates nodes if there are omitted element blocks.
  if (interFace.match_node_ids()) {
    build_reverse_node_map(output_region, part_mesh, global_node_map, local_node_map);
  }
  else if (interFace.match_node_xyz()) {
    match_node_xyz(part_mesh, interFace.tolerance(), global_node_map, local_node_map);
  }
  else if (interFace.match_nodeset_nodes()) {
    match_nodeset_nodes(part_mesh, interFace.tolerance(), global_node_map, local_node_map,
                        interFace);
  }
  else {
    // Eliminate all nodes that were only connected to the omitted element blocks (if any).
    bool fill_global = true;
    eliminate_omitted_nodes(part_mesh, global_node_map, local_node_map, fill_global);
  }

  node_count    = global_node_map.size();
  size_t merged = local_node_map.size() - global_node_map.size();
  if (merged > 0) {
    fmt::print("*** {} Nodes were merged/omitted.\n", fmt::group_digits(merged));
  }

// Verify nodemap...
#ifndef NDEBUG
  std::vector<int> glob(node_count);
  for (auto id : local_node_map) {
    if (id >= 0) {
      glob[id] = 1;
    }
  }
  for (int i : glob) {
    SMART_ASSERT(i == 1);
  }
#endif

  // Transfer some common data...
  output_region.property_add(part_mesh[0]->get_property("title"));

  // Define a node block...
  std::string block_name        = "nodeblock_1";
  int         spatial_dimension = part_mesh[0]->get_property("spatial_dimension").get_int();
  auto       *block =
      new Ioss::NodeBlock(output_region.get_database(), block_name, node_count, spatial_dimension);
  block->property_add(Ioss::Property("id", 1));

  output_region.add(block);

  const std::string &eb_combines = interFace.elementblock_combines();
  if (!eb_combines.empty()) {
    process_specified_combines(part_mesh, eb_combines, Ioss::ELEMENTBLOCK);
  }

  const std::string ss_combines = interFace.sideset_combines();
  if (!ss_combines.empty()) {
    process_specified_combines(part_mesh, ss_combines, Ioss::SIDESET);
  }

  const std::string ns_combines = interFace.nodeset_combines();
  if (!ns_combines.empty()) {
    process_specified_combines(part_mesh, ns_combines, Ioss::NODESET);
  }

  // Add element blocks, nodesets, sidesets
  for (size_t p = 0; p < part_count; p++) {
    transfer_elementblock(*part_mesh[p], output_region, interFace.create_assemblies(),
                          interFace.combine_element_blocks(), false);
    if (interFace.convert_nodes_to_nodesets(p + 1)) {
      create_nodal_nodeset(*part_mesh[p], output_region, false);
    }
    if (!interFace.omit_nodesets()) {
      transfer_nodesets(*part_mesh[p], output_region, interFace.combine_nodesets(), false);
    }
    if (!interFace.omit_sidesets()) {
      transfer_sidesets(*part_mesh[p], output_region, interFace.combine_sidesets(), false);
    }
    if (!interFace.omit_assemblies()) {
      transfer_assembly(*part_mesh[p], output_region, false);
    }
  }

  if (merged > 0 && interFace.combine_nodesets() &&
      (interFace.match_node_xyz() || interFace.match_nodeset_nodes() ||
       interFace.match_node_ids())) {
    // Get the nodelist for each combined nodeset and see if contains duplicate nodes...
    check_for_duplicate_nodeset_nodes(output_region, local_node_map);
  }

  // This is the map from local element position to global element
  // position (0-based). If there are no element block omissions or
  // combinations, then the map is simply [0..number_elements). If
  // there are omissions, Then local_element_map[j] will be -1 for an
  // omitted element.  If there are element block combinations, then
  // the map will not be 0..number_elements. If not omitted,
  // local_element_map[part_offset+j] gives the global position of
  // local element j in the current part.
  //
  // This needs to be constructed after the element blocks are
  // combined.
  std::vector<INT> local_element_map =
      build_local_element_map<INT>(part_mesh, output_region, output_input_map);

  if (!interFace.information_record_parts().empty()) {
    const std::vector<int> &info_parts = interFace.information_record_parts();
    if (info_parts[0] == 0) {
      // Transfer info records from all parts...
      for (const auto &pm : part_mesh) {
        const StringVector &info = pm->get_information_records();
        output_region.add_information_records(info);
      }
    }
    else {
      for (int info_part : info_parts) {
        const StringVector &info = part_mesh[info_part - 1]->get_information_records();
        output_region.add_information_records(info);
      }
    }
  }

  output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

  output_region.begin_mode(Ioss::STATE_MODEL);

  output_nodeblock(output_region, part_mesh, local_node_map, global_node_map);
  output_elementblock(output_region, part_mesh, local_node_map, local_element_map,
                      interFace.ignore_element_ids() || interFace.combine_element_blocks());
  output_nodal_nodeset(output_region, part_mesh, interFace, local_node_map);

  if (!interFace.omit_nodesets()) {
    output_nodeset(output_region, local_node_map, merged > 0);
  }
  if (!interFace.omit_sidesets()) {
    output_sideset(output_region, local_element_map);
  }

  output_region.end_mode(Ioss::STATE_MODEL);

  // ####################TRANSIENT DATA SECTION###########################
  // ***********************************************************************
  // 9. Get Variable Information and names

  if (debug_level & 1) {
    fmt::print(stderr, "{}", time_stamp(tsFormat));
  }

  output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

  bool error = false;
  error |= define_global_fields(output_region, part_mesh, interFace.global_var_names());

  error |= define_nodal_fields(output_region, part_mesh, interFace.node_var_names(), interFace);
  define_nodal_nodeset_fields(output_region, part_mesh, interFace.node_var_names(), interFace);

  error |= define_element_fields(output_region, interFace.elem_var_names());

  if (!interFace.omit_nodesets()) {
    error |= define_nodeset_fields(output_region, interFace.nodeset_var_names());
  }
  if (!interFace.omit_sidesets()) {
    error |= define_sideset_fields(output_region, interFace.sideset_var_names());
  }

  output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

  if (error) {
    fmt::print(stderr,
               "ERROR: Specified field(s) (see above) were not found. Fix input and rerun.\n\n");
    exit(EXIT_FAILURE);
  }

  // Get database times...
  // Different parts can have either no times or the times must match
  //! \todo If ts_min, ts_max, ts_step is specified, then only check steps in that range...
  std::vector<double> global_times;

  for (size_t p = 0; p < part_count; p++) {
    int                 nts = part_mesh[p]->get_property("state_count").get_int();
    std::vector<double> times(nts);
    // If multiple databases have timesteps, they must all match.
    for (int i = 0; i < nts; i++) {
      times[i] = part_mesh[p]->get_state_time(i + 1);
    }

    if (nts > 0) {
      if (global_times.empty()) {
        std::copy(times.begin(), times.end(), std::back_inserter(global_times));
      }
      else {
        if (global_times.size() != times.size()) {
          fmt::print(stderr, "ERROR: Time step sizes must match. (Global = {}, Part {} = {}",
                     global_times.size(), p, times.size());
          exit(EXIT_FAILURE);
        }

        for (size_t i = 0; i < global_times.size(); i++) {
          if (!approx_equal(global_times[i], times[i], global_times[0])) {
            fmt::print(stderr,
                       "ERROR: Time step {} in part {} does not match time steps in previous "
                       "part(s): previous: "
                       "{}, current: {}\n",
                       i + 1, p + 1, global_times[i], times[i]);
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  }

  output_region.begin_mode(Ioss::STATE_TRANSIENT);
  fmt::print("\n");
  int ts_min    = interFace.step_min();
  int ts_max    = interFace.step_max();
  int ts_step   = interFace.step_interval();
  int num_steps = static_cast<int>(global_times.size());

  if (ts_min < 0) {
    ts_min = num_steps + 1 + ts_min;
  }
  if (ts_max < 0) {
    ts_max = num_steps + 1 + ts_max;
  }
  ts_max = ts_max < num_steps ? ts_max : num_steps;

  double ts_begin = Ioss::Utils::timer();
  int    steps    = 0;
  int    nsteps   = (ts_max - ts_min + 1) / ts_step;
  for (int step = ts_min - 1; step < ts_max; step += ts_step) {
    int ostep = output_region.add_state(global_times[step]);
    output_region.begin_state(ostep);
    output_transient_state(output_region, part_mesh, global_times[step], local_node_map, interFace,
                           merged);
    fmt::print("\rWrote step {:4}/{:4}, time {:8.4e}", step + 1, nsteps, global_times[step]);
    output_region.end_state(ostep);
    steps++;
  }
  double end = Ioss::Utils::timer();
  fmt::print("\n");
  output_region.end_mode(Ioss::STATE_TRANSIENT);

  /*************************************************************************/
  // EXIT program
  if (debug_level & 1) {
    fmt::print(stderr, "{}", time_stamp(tsFormat));
  }
  output_region.output_summary(std::cout);
  fmt::print("******* END *******\n");
  fmt::print(stderr, "\nTotal Execution Time     = {:.5} seconds.\n", end - begin);
  if (steps > 0) {
    fmt::print(stderr, "\tMesh = {:.5} seconds; Timesteps = {:.5} seconds / step.\n\n",
               (ts_begin - begin), (end - ts_begin) / (double)(steps));
  }
  return end - begin;
}

namespace {
  bool entity_is_omitted(Ioss::GroupingEntity *block)
  {
    bool omitted = block->get_optional_property("omitted", 0) == 1;
    return omitted;
  }

  void transfer_elementblock(const Ioss::Region &region, Ioss::Region &output_region,
                             bool create_assemblies, bool combine_similar, bool /* debug */)
  {
    static int         used_blocks = 0;
    const std::string &prefix      = region.name();

    Ioss::Assembly *assem = nullptr;

    const Ioss::ElementBlockContainer &ebs = region.get_element_blocks();
    for (const auto &eb : ebs) {
      if (!entity_is_omitted(eb)) {
        std::string name = eb->name();
        name             = eb->get_optional_property("ejoin_combine_into", name);
        auto *oeb        = output_region.get_element_block(name);
        if (oeb != nullptr) {
          if (combine_similar || eb->property_exists("ejoin_combine_into")) {
            if (oeb->topology() != eb->topology()) {
              fmt::print(
                  stderr,
                  "ERROR: The topology ('{}') for element block '{}' does not match\n       the "
                  "topology ('{}') for element block '{}'.\n       They cannot be combined.\n\n",
                  oeb->topology()->name(), oeb->name(), eb->topology()->name(), eb->name());
              exit(EXIT_FAILURE);
            }
            // Combine element blocks with similar names...
            output_input_map[oeb].emplace_back(eb, oeb->entity_count());
            size_t count = eb->entity_count();
            add_to_entity_count(oeb, count);
            continue;
          }
          else {
            name = prefix + "_" + eb->name();
            if (output_region.get_element_block(name) != nullptr) {
              fmt::print(stderr, "ERROR: Duplicate element blocks named '{}'\n", name);
              exit(EXIT_FAILURE);
            }
          }
        }
        // This is a new element block at this point...
        eb->property_add(Ioss::Property("name_in_output", name));

        std::string type     = eb->topology()->name();
        size_t      num_elem = eb->entity_count();

        if (num_elem > 0) {
          auto *ebn = new Ioss::ElementBlock(output_region.get_database(), name, type, num_elem);
          output_input_map[ebn].emplace_back(eb, 0);
          ebn->property_add(Ioss::Property("original_block_order", used_blocks++));
          output_region.add(ebn);
          transfer_fields(eb, ebn, Ioss::Field::ATTRIBUTE);

          set_id(eb, ebn);

          if (eb->property_exists("original_topology_type")) {
            std::string oes = eb->get_property("original_topology_type").get_string();

            // Set the new property
            ebn->property_add(Ioss::Property("original_topology_type", oes));
          }

          if (create_assemblies) {
            if (assem == nullptr) {
              assem = new Ioss::Assembly(output_region.get_database(), region.name());
              output_region.add(assem);
              assem->property_add(Ioss::Property("name_in_output", region.name()));
            }
            assem->add(ebn);
          }
        }
      }
    }
  }

  void transfer_assembly(const Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    // All assemblies on the input parts will be transferred to the output mesh
    // Possibly renamed if a name conflict
    // Also need to update names of the entities in the assembly since they were possibly
    // renamed on the output...
    // TODO: handle renamed nested assemblies...
    const std::string &prefix = region.name();

    const Ioss::AssemblyContainer &assems = region.get_assemblies();
    for (const auto &as : assems) {
      if (!entity_is_omitted(as)) {
        std::string name = as->name();
        if (output_region.get_assembly(name) != nullptr) {
          name = prefix + "_" + as->name();
          if (output_region.get_assembly(name) != nullptr) {
            fmt::print(stderr, "ERROR: Duplicate assemblies named '{}'\n", name);
            exit(EXIT_FAILURE);
          }
        }
        as->property_add(Ioss::Property("name_in_output", name));

        if (debug) {
          fmt::print(stderr, "{}, ", name);
        }
        size_t num_members = as->entity_count();
        if (num_members > 0) {
          auto  member_type = as->get_member_type();
          auto *asn         = new Ioss::Assembly(output_region.get_database(), name);
          output_region.add(asn);

          const auto &members = as->get_members();
          for (const auto &member : members) {
            auto output_name = member->get_property("name_in_output").get_string();

            auto *entity = output_region.get_entity(output_name, member_type);
            SMART_ASSERT(entity != nullptr);
            asn->add(entity);
          }
        }
      }
    }
  }

  void transfer_sideblocks(const Ioss::Region &region, Ioss::Region &output_region,
                           bool combine_similar)
  {
    const std::string            &prefix = region.name();
    const Ioss::SideSetContainer &osss   = output_region.get_sidesets();

    for (const auto &oss : osss) {
      const auto &itr = output_input_map.find(oss);
      SMART_ASSERT(itr != output_input_map.end());
      const auto &[key, oss_inputs] = *itr;
      for (const auto &[gss, offset] : oss_inputs) {
        if (gss != nullptr) {
          auto *iss = dynamic_cast<const Ioss::SideSet *>(gss);
          Ioss::Utils::check_dynamic_cast(iss);
          if (*(iss->contained_in()) == region) {
            const Ioss::SideBlockContainer &sbs = iss->get_side_blocks();

            for (const auto &sb : sbs) {
              std::string sbname = sb->name();
              auto       *osb    = oss->get_side_block(sbname);
              if (osb != nullptr) {
                if (combine_similar) {
                  // Combine side blocks with similar names...
                  output_input_map[osb].emplace_back(sb, osb->entity_count());
                  size_t count = sb->entity_count();
                  add_to_entity_count(osb, count);
                  continue;
                }
                else {
                  sbname = prefix + "_" + sb->name();
                  if (oss->get_side_block(sbname) != nullptr) {
                    fmt::print(stderr, "ERROR: Duplicate sideset sideblocks named '{}'\n", sbname);
                    exit(EXIT_FAILURE);
                  }
                }
              }
              // This is a new sideblock at this point...
              sb->property_add(Ioss::Property("name_in_output", sbname));
              std::string sbtype   = sb->topology()->name();
              std::string partype  = sb->parent_element_topology()->name();
              size_t      num_side = sb->entity_count();
              auto       *block = new Ioss::SideBlock(output_region.get_database(), sbname, sbtype,
                                                      partype, num_side);
              output_input_map[block].emplace_back(sb, 0);
              oss->add(block);
            }
          }
        }
      }
    }
  }

  void transfer_sidesets(const Ioss::Region &region, Ioss::Region &output_region,
                         bool combine_similar, bool /* debug */)
  {
    const std::string &prefix = region.name();

    const Ioss::SideSetContainer &sss = region.get_sidesets();
    for (const auto &ss : sss) {
      if (!entity_is_omitted(ss)) {
        std::string name = ss->name();
        name             = ss->get_optional_property("ejoin_combine_into", name);
        auto *oss        = output_region.get_sideset(name);
        if (oss != nullptr) {
          if (combine_similar || ss->property_exists("ejoin_combine_into")) {
            // Combine sidesets with similar names...
            output_input_map[oss].emplace_back(ss, oss->entity_count());
            size_t count = ss->entity_count();
            add_to_entity_count(oss, count);
            continue;
          }
          else {
            name = prefix + "_" + ss->name();
            if (output_region.get_sideset(name) != nullptr) {
              fmt::print(stderr, "ERROR: Duplicate side sets named '{}'\n", name);
              exit(EXIT_FAILURE);
            }
          }
        }
        // This is a new sideset at this point...
        ss->property_add(Ioss::Property("name_in_output", name));
        auto *surf = new Ioss::SideSet(output_region.get_database(), name);
        output_input_map[surf].emplace_back(ss, 0);
        set_id(ss, surf);
        output_region.add(surf);
      }
    }

    // Now deal with the sideblocks...
    transfer_sideblocks(region, output_region, combine_similar);
  }

  // Create a nodeset on the output region consisting of all the nodes
  // in the input region.
  void create_nodal_nodeset(const Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const std::string &prefix = region.name();

    std::string name = prefix + "_nodes";
    if (output_region.get_nodeset(name) != nullptr) {
      fmt::print(stderr, "ERROR: Duplicate node sets named '{}'\n", name);
      exit(EXIT_FAILURE);
    }
    if (debug) {
      fmt::print(stderr, "{}, ", name);
    }
    size_t count = region.get_property("node_count").get_int();
    auto  *ns    = new Ioss::NodeSet(output_region.get_database(), name, count);
    output_region.add(ns);
    output_input_map[ns].emplace_back(nullptr, 0);
  }

  // Output the bulk data for a nodeset on the output region
  // consisting of all the nodes in the input region.
  template <typename INT>
  void output_nodal_nodeset(Ioss::Region &output_region, const RegionVector &part_mesh,
                            const SystemInterface  &interFace,
                            const std::vector<INT> &local_node_map)
  {
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      if (interFace.convert_nodes_to_nodesets(p + 1)) {
        std::string    name = part_mesh[p]->name() + "_nodes";
        Ioss::NodeSet *ons  = output_region.get_nodeset(name);
        SMART_ASSERT(ons != nullptr)(name);

        Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
        SMART_ASSERT(nb != nullptr);
        size_t node_offset = part_mesh[p]->get_property("node_offset").get_int();

        std::vector<INT> nodelist;
        nb->get_field_data("ids", nodelist);
        // This needs to make sure that the nodelist comes back as local id (1..numnodes)
        for (auto &node : nodelist) {
          size_t loc_node = part_mesh[p]->node_global_to_local(node, true) - 1;
          auto   gpos     = local_node_map[node_offset + loc_node];
          if (gpos >= 0) {
            node = gpos + 1;
          }
        }
        ons->put_field_data("ids_raw", nodelist);

        // Output distribution factors -- set all to 1.0
        std::vector<double> factors(nodelist.size(), 1.0);
        ons->put_field_data("distribution_factors", factors);
      }
    }
  }

  void define_nodal_nodeset_fields(const Ioss::Region &output_region, const RegionVector &part_mesh,
                                   const StringIdVector  &variable_list,
                                   const SystemInterface &interFace)
  {
    // This routine does not check that all variables in `variable_list` have been
    // found since the checking has already been done in define_nodal_fields.
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return;
    }

    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {

      if (interFace.convert_nodes_to_nodesets(p + 1)) {
        // Find nodeset in output region corresponding to the nodes in this part...
        std::string    name = part_mesh[p]->name() + "_nodes";
        Ioss::NodeSet *ons  = output_region.get_nodeset(name);
        SMART_ASSERT(ons != nullptr)(name);

        Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
        SMART_ASSERT(nb != nullptr);

        SMART_ASSERT(part_mesh[p]->get_property("node_count").get_int() == nb->entity_count());

        Ioss::NameList fields = nb->field_describe(Ioss::Field::TRANSIENT);
        for (const auto &field_name : fields) {
          if (valid_variable(field_name, 0, variable_list)) {
            Ioss::Field field = nb->get_field(field_name);
            ons->field_add(std::move(field));
          }
        }
      }
    }
  }

  void transfer_nodesets(const Ioss::Region &region, Ioss::Region &output_region,
                         bool combine_similar, bool /*debug*/)
  {
    const std::string &prefix = region.name();

    const Ioss::NodeSetContainer &nss = region.get_nodesets();
    for (const auto &ns : nss) {
      if (!entity_is_omitted(ns)) {
        std::string name = ns->name();
        name             = ns->get_optional_property("ejoin_combine_into", name);
        auto *ons        = output_region.get_nodeset(name);
        if (ons != nullptr) {
          if (combine_similar || ns->property_exists("ejoin_combine_into")) {
            // Combine nodesets with similar names...
            output_input_map[ons].emplace_back(ns, ons->entity_count());
            size_t count = ns->entity_count();
            add_to_entity_count(ons, count);
            continue;
          }
          else {
            name = prefix + "_" + ns->name();
            if (output_region.get_nodeset(name) != nullptr) {
              fmt::print(stderr, "ERROR: Duplicate node sets named '{}'\n", name);
              exit(EXIT_FAILURE);
            }
          }
        }
        // This is a new output nodeset at this point...
        size_t count    = ns->entity_count();
        auto  *node_set = new Ioss::NodeSet(output_region.get_database(), name, count);
        output_input_map[node_set].emplace_back(ns, 0);
        output_region.add(node_set);
        set_id(ns, node_set);
      }
    }
  }

  template <typename T, typename INT>
  void map_element_vars(size_t loffset, size_t goffset, size_t entity_count, std::vector<T> &values,
                        std::vector<T> &global_values, INT *part_loc_elem_to_global)
  {
    // copy values to master element value information
    T *local_values = Data(values);
    for (size_t j = 0; j < entity_count; j++) {
      size_t global_block_pos         = part_loc_elem_to_global[(j + loffset)] - goffset;
      global_values[global_block_pos] = local_values[j];
    }
  }

  template <typename T>
  void map_sideset_vars(size_t loffset, size_t entity_count, std::vector<T> &values,
                        std::vector<T> &global_values)
  {
    // copy values to master sideset value information
    T *local_values = Data(values);
    for (size_t j = 0; j < entity_count; j++) {
      global_values[j + loffset] = local_values[j];
    }
  }

  template <typename T, typename U>
  void map_nodeset_vars(U & /*unused*/, int /*unused*/, int /*unused*/, std::vector<T> & /*unused*/,
                        std::vector<T> & /*unused*/)
  {
    SMART_ASSERT(1 == 0 && "Internal Error!");
  }
} // namespace

namespace {
  template <typename INT>
  void output_nodeblock(Ioss::Region &output_region, const RegionVector &part_mesh,
                        const std::vector<INT> &local_node_map, std::vector<INT> &global_node_map)
  {
    Ioss::NodeBlock *onb = output_region.get_node_blocks()[0];
    SMART_ASSERT(onb != nullptr);

    onb->put_field_data("ids", global_node_map);

    int spatial_dimension = output_region.get_property("spatial_dimension").get_int();
    std::vector<double> coord(global_node_map.size() * spatial_dimension);
    for (const auto &pm : part_mesh) {
      Ioss::NodeBlock *nb = pm->get_node_blocks()[0];
      SMART_ASSERT(nb != nullptr);
      std::vector<double> coordinates;
      nb->get_field_data("mesh_model_coordinates", coordinates);
      size_t node_count = nb->entity_count();
      size_t offset     = pm->get_property("node_offset").get_int();
      for (size_t i = 0; i < node_count; i++) {
        auto glob_pos = local_node_map[i + offset];
        if (glob_pos >= 0) {
          coord[glob_pos * spatial_dimension + 0] = coordinates[i * spatial_dimension + 0];
          coord[glob_pos * spatial_dimension + 1] = coordinates[i * spatial_dimension + 1];
          coord[glob_pos * spatial_dimension + 2] = coordinates[i * spatial_dimension + 2];
        }
      }
    }
    onb->put_field_data("mesh_model_coordinates", coord);
  }

  template <typename INT>
  void output_elementblock(Ioss::Region &output_region, const RegionVector &part_mesh,
                           const std::vector<INT> &local_node_map,
                           const std::vector<INT> &local_element_map, bool ignore_element_ids)
  {

    const Ioss::ElementBlockContainer &ebs = output_region.get_element_blocks();

    // Ids...
    size_t           element_count = output_region.get_property("element_count").get_int();
    std::vector<INT> ids(element_count);

    if (ignore_element_ids) {
      // Just generate 1..numel ids (much faster for large models)
      std::iota(ids.begin(), ids.end(), 1);
    }
    else {
      // Try to maintain the original element ids if possible...
      generate_element_ids(part_mesh, local_element_map, ids);
    }
    size_t element_offset = 0;
    for (const auto &eb : ebs) {
      eb->put_field_data("ids", &ids[element_offset], ids.size() * sizeof(INT));
      element_offset += eb->entity_count();
    }

    SMART_ASSERT(element_offset == element_count);

    // Connectivity...
    for (const auto &oeb : ebs) {
      const auto &itr = output_input_map.find(oeb);
      SMART_ASSERT(itr != output_input_map.end());
      const auto &[key, oeb_inputs] = *itr;
      if (!oeb_inputs.empty()) {
        int64_t          count = oeb->entity_count();
        int64_t          nnpe  = oeb->topology()->number_nodes();
        std::vector<INT> connectivity(count * nnpe);
        for (const auto &[ieb, offset] : oeb_inputs) {
          if (ieb != nullptr) {
            ieb->get_field_data("connectivity_raw", &connectivity[offset * nnpe], -1);

            auto *input_region = dynamic_cast<const Ioss::Region *>(ieb->contained_in());
            Ioss::Utils::check_dynamic_cast(input_region);
            size_t node_offset = input_region->get_property("node_offset").get_int();
            for (int64_t i = 0; i < ieb->entity_count() * nnpe; i++) {
              // connectivity is in part-local node ids [1..num_node]
              // loc_node = the position of node in the local [0..num_node)
              // local_node_map[node_offset+loc_node] gives the position of this node in the
              // global list
              size_t loc_node = connectivity[offset * nnpe + i] - 1;
              SMART_ASSERT(node_offset + loc_node < local_node_map.size());
              auto gpos = local_node_map[node_offset + loc_node];
              if (gpos >= 0) {
                connectivity[offset * nnpe + i] = gpos + 1;
              }
            }
          }
        }
        oeb->put_field_data("connectivity_raw", connectivity);
      }
    }
  }

  template <typename INT, typename T> size_t unique(std::vector<std::pair<INT, T>> &out)
  {
    if (out.empty()) {
      return 0;
    }
    size_t            i    = 1;
    size_t            pos  = 1;
    std::pair<INT, T> oldv = out[0];
    for (; i < out.size(); ++i) {
      std::pair<INT, T> newv = out[i];
      out[pos]               = newv;
      pos += (newv.first != oldv.first);
      oldv = newv;
    }
    out.resize(pos);
    out.shrink_to_fit();
    return pos;
  }

  template <typename INT>
  void output_nodeset(Ioss::Region &output_region, const std::vector<INT> &local_node_map,
                      bool nodes_consolidated)
  {
    const auto &output_nodesets = output_region.get_nodesets();
    for (const auto &ons : output_nodesets) {
      const auto &itr = output_input_map.find(ons);
      SMART_ASSERT(itr != output_input_map.end());
      const auto &[key, ons_inputs] = *itr;
      if (!ons_inputs.empty()) {
        int64_t count = ons->entity_count();

        // The size of the input nodeset nodelists may be more than the size of the output nodeset
        // nodelist due to node consolidation...
        if (nodes_consolidated) {
          count = 0;
          for (const auto &[ins, offset] : ons_inputs) {
            if (ins != nullptr) {
              count += ins->entity_count();
            }
          }
        }

        std::vector<INT>    nodelist(count);
        std::vector<double> df(count);
        bool                found_one = false;
        for (const auto &[ins, offset] : ons_inputs) {
          if (ins != nullptr) {
            found_one = true;
            ins->get_field_data("ids", &nodelist[offset], -1);

            auto *input_region = dynamic_cast<const Ioss::Region *>(ins->contained_in());
            Ioss::Utils::check_dynamic_cast(input_region);
            size_t node_offset = input_region->get_property("node_offset").get_int();
            for (int64_t i = 0; i < ins->entity_count(); i++) {
              size_t loc_node = input_region->node_global_to_local(nodelist[offset + i], true) - 1;
              auto   gpos     = local_node_map[node_offset + loc_node];
              if (gpos >= 0) {
                nodelist[offset + i] = gpos + 1;
              }
            }
            ins->get_field_data("distribution_factors", &df[offset], -1);
          }
        }
        if (found_one) {
          if (count != ons->entity_count()) {
            const auto &ns_itr = nodeset_in_out_map.find(ons);
            SMART_ASSERT(ns_itr != nodeset_in_out_map.end());
            const auto &[ns_key, map] = *ns_itr;
            SMART_ASSERT(ns_key == ons);
            SMART_ASSERT(map.size() == (size_t)ons->entity_count());
            for (int64_t i = 0; i < ons->entity_count(); i++) {
              nodelist[i] = nodelist[map[i]];
              df[i]       = df[map[i]];
            }
          }
          ons->put_field_data("ids_raw", nodelist);
          ons->put_field_data("distribution_factors", df);
        }
      }
    }
  }

  template <typename INT>
  void output_sideblock(Ioss::SideSet *oss, const std::vector<INT> &local_element_map)
  {
    // Get output sideblocks in the output sideset `oss`
    const Ioss::SideBlockContainer &osbs = oss->get_side_blocks();
    for (const auto &osb : osbs) {
      const auto &itr = output_input_map.find(osb);
      SMART_ASSERT(itr != output_input_map.end());
      const auto &[key, osb_inputs] = *itr;
      if (!osb_inputs.empty()) {
        int64_t          count = osb->entity_count();
        std::vector<INT> elem_side_list(count * 2);

        for (const auto &[isb, offset] : osb_inputs) {
          if (isb != nullptr) {
            isb->get_field_data("element_side_raw", &elem_side_list[offset * 2], -1);
            auto *input_region = isb->contained_in()->contained_in();
            SMART_ASSERT(input_region != nullptr);
            size_t element_offset = input_region->get_property("element_offset").get_int();

            // The 'elem_side_list' contains
            // (local_element_position,side_ordinal) pairs. The
            // 'local_element_position' is 1-based offset in the
            // current part.  Need to map to its location in the
            // output region...
            for (int64_t i = 0; i < isb->entity_count() * 2;
                 i += 2) { // just modify the elem part of the pair...
              size_t local_position = elem_side_list[2 * offset + i] - 1;
              auto   gpos           = local_element_map[element_offset + local_position];
              SMART_ASSERT(gpos >= 0)(gpos)(i); // Inactive elements should be filtered by Ioss
              elem_side_list[2 * offset + i] = gpos + 1;
            }
          }
        }
        osb->put_field_data("element_side_raw", elem_side_list);
      }
    }
  }

  template <typename INT>
  void output_sideset(Ioss::Region &output_region, const std::vector<INT> &local_element_map)
  {
    const auto &output_sidesets = output_region.get_sidesets();
    for (const auto &oss : output_sidesets) {
      output_sideblock(oss, local_element_map);
    }
  }

  void output_global_fields(Ioss::Region &output_region, const RegionVector &part_mesh)
  {
    for (const auto &pm : part_mesh) {
      Ioss::NameList fields = pm->field_describe(Ioss::Field::REDUCTION);
      for (const auto &field : fields) {
        std::vector<double> data;
        pm->get_field_data(field, data);
        output_region.put_field_data(field, data);
      }
    }
  }

  template <typename INT>
  void output_nodal_fields(Ioss::Region &output_region, const RegionVector &part_mesh,
                           const std::vector<INT> &local_node_map, SystemInterface &interFace)
  {
    size_t part_count = part_mesh.size();

    Ioss::NodeBlock *onb = output_region.get_node_blocks()[0];
    SMART_ASSERT(onb != nullptr);
    size_t node_count = onb->entity_count();

    Ioss::NameList fields = onb->field_describe(Ioss::Field::TRANSIENT);
    for (const auto &field : fields) {
      size_t              comp_count = onb->get_field(field).raw_storage()->component_count();
      std::vector<double> data(node_count * comp_count);
      for (size_t p = 0; p < part_count; p++) {
        if (!interFace.convert_nodes_to_nodesets(p + 1)) {
          size_t           offset = part_mesh[p]->get_property("node_offset").get_int();
          Ioss::NodeBlock *nb     = part_mesh[p]->get_node_blocks()[0];
          SMART_ASSERT(nb != nullptr);
          if (nb->field_exists(field)) {
            SMART_ASSERT((int)comp_count == nb->get_field(field).raw_storage()->component_count());
            std::vector<double> loc_data;
            nb->get_field_data(field, loc_data);
            size_t nc = nb->entity_count();
            SMART_ASSERT(loc_data.size() == nc * comp_count);
            for (size_t i = 0; i < nc; i++) {
              auto glob_pos = local_node_map[offset + i];
              if (glob_pos >= 0) {
                for (size_t j = 0; j < comp_count; j++) {
                  data[glob_pos * comp_count + j] = loc_data[i * comp_count + j];
                }
              }
            }
          }
        }
      }
      onb->put_field_data(field, data);
    }
  }

  void output_nodal_nodeset_fields(Ioss::Region &output_region, const RegionVector &part_mesh,
                                   SystemInterface &interFace)
  {
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      if (interFace.convert_nodes_to_nodesets(p + 1)) {

        // Find nodeset in output region corresponding to the nodes in this part...
        std::string    name = part_mesh[p]->name() + "_nodes";
        Ioss::NodeSet *ons  = output_region.get_nodeset(name);
        SMART_ASSERT(ons != nullptr)(name);

        Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
        SMART_ASSERT(nb != nullptr);

        SMART_ASSERT(part_mesh[p]->get_property("node_count").get_int() == nb->entity_count());

        // NOTE: The node order in the output nodeset 'ons' was
        // defined as the same node order in the input nodeblock 'nb',
        // so we shouldn't have to do any reordering of the data at
        // this time--just read then write.
        Ioss::NameList      fields = ons->field_describe(Ioss::Field::TRANSIENT);
        std::vector<double> data;
        for (const auto &field : fields) {
          nb->get_field_data(field, data);
          ons->put_field_data(field, data);
        }
      }
    }
  }

  void output_entity_fields(Ioss::GroupingEntity *out_entity)
  {
    const auto &itr = output_input_map.find(out_entity);
    if (itr != output_input_map.end()) {
      const auto &[key, out_entity_inputs] = *itr;
      if (!out_entity_inputs.empty()) {
        Ioss::NameList fields = out_entity->field_describe(Ioss::Field::TRANSIENT);
        int64_t        count  = out_entity->entity_count();

        for (const auto &field : fields) {
          bool    found_field = false;
          int64_t comp_count =
              out_entity->get_field(field).get_component_count(Ioss::Field::InOut::OUTPUT);
          std::vector<double> field_data(comp_count * count);
          for (const auto &[ieb, offset] : out_entity_inputs) {
            if (ieb != nullptr && ieb->field_exists(field)) {
              found_field = true;
              ieb->get_field_data(field, &field_data[comp_count * offset], -1);
            }
          }
          if (found_field) {
            out_entity->put_field_data(field, field_data);
          }
        }
      }
    }
  }

  void output_element_fields(Ioss::Region &output_region)
  {
    const auto &output_element_blocks = output_region.get_element_blocks();
    for (const auto &oeb : output_element_blocks) {
      output_entity_fields(oeb);
    }
  }

  void output_nodeset_fields(Ioss::Region &output_region, bool nodes_consolidated)
  {
    // NOTE: The handling of merged nodes is very inefficient currently since it is done once per
    // timestep...
    //       See if can store the map somewhere and use it.  It is initially built in
    //       output_nodeset...
    const auto &output_nodesets = output_region.get_nodesets();
    if (nodes_consolidated) {
      // The size of the input nodeset nodelists may be more than the size of the output nodeset
      // nodelist due to node consolidation...
      for (const auto &ons : output_nodesets) {
        const auto &itr = output_input_map.find(ons);
        SMART_ASSERT(itr != output_input_map.end());
        const auto &[key, ons_inputs] = *itr;
        int64_t count                 = 0;
        for (const auto &[ins, offset] : ons_inputs) {
          if (ins != nullptr) {
            count += ins->entity_count();
          }
        }

        if (count == 0) {
          continue;
        }
        if (count == ons->entity_count()) {
          output_entity_fields(ons);
        }
        else {
          // There is at least one duplicated node that is removed.  Need to map nodeset fields to
          // account for deleted node(s)
          Ioss::NameList fields = ons->field_describe(Ioss::Field::TRANSIENT);
          if (fields.empty()) {
            continue;
          }

          // Get the mapping of the input nodelist node position to the output position...
          // If this is the `nodal_nodeset`, then it will not be found in `nodeset_in_out_map`
          const auto &ns_itr = nodeset_in_out_map.find(ons);
          if (ns_itr != nodeset_in_out_map.end()) {
            const auto &[ns_key, map] = *ns_itr;
            SMART_ASSERT(ns_key == ons);
            SMART_ASSERT(map.size() == (size_t)ons->entity_count());

            // Now get each field, map to correct output position and output...
            for (const auto &field_name : fields) {
              size_t comp_count = ons->get_field(field_name).raw_storage()->component_count();
              std::vector<double> field_data(count * comp_count);

              for (const auto &[ins, offset] : ons_inputs) {
                if (ins != nullptr && ins->field_exists(field_name)) {
                  ins->get_field_data(field_name, &field_data[comp_count * offset], -1);
                }
              }
              for (int64_t i = 0; i < ons->entity_count(); i++) {
                for (size_t j = 0; j < comp_count; j++) {
                  field_data[comp_count * i + j] = field_data[comp_count * map[i] + j];
                }
              }
              ons->put_field_data(field_name, field_data);
            }
          }
        }
      }
    }
    else {
      for (const auto &ons : output_nodesets) {
        output_entity_fields(ons);
      }
    }
  }

  void output_sideblock_fields(Ioss::SideSet *output_sideset)
  {
    const Ioss::SideBlockContainer &output_sideblocks = output_sideset->get_side_blocks();
    for (const auto &osb : output_sideblocks) {
      output_entity_fields(osb);
    }
  }

  void output_sideset_fields(Ioss::Region &output_region)
  {
    const Ioss::SideSetContainer &output_sidesets = output_region.get_sidesets();
    for (const auto &oss : output_sidesets) {
      const auto &itr = output_input_map.find(oss);
      SMART_ASSERT(itr != output_input_map.end());
      const auto &[key, oss_inputs] = *itr;
      if (!oss_inputs.empty()) {
        output_sideblock_fields(oss);
      }
    }
  }

  template <typename INT>
  void output_transient_state(Ioss::Region &output_region, const RegionVector &part_mesh,
                              double time, const std::vector<INT> &local_node_map,
                              SystemInterface &interFace, bool merged)
  {
    // Determine which state on each input mesh corresponds to 'time'
    std::vector<int> steps(part_mesh.size());
    for (size_t p = 0; p < part_mesh.size(); p++) {
      double min_delta = 1.0e39;
      size_t min_step  = 0;

      size_t nts = part_mesh[p]->get_property("state_count").get_int();
      steps[p]   = 0;
      for (size_t i = 0; i < nts; i++) {
        double delta = std::fabs(part_mesh[p]->get_state_time(i + 1) - time);
        if (delta < min_delta) {
          min_delta = delta;
          min_step  = i;
          if (delta == 0.0) {
            break;
          }
        }
        else {
          // Delta is increasing; times are moving apart...
          // (Assumes monotonically increasing time values...)
          break;
        }
      }
      if (nts > 0) {
        steps[p] = min_step + 1;
        part_mesh[p]->begin_state(steps[p]);
      }
    }

    output_global_fields(output_region, part_mesh);
    output_nodal_fields(output_region, part_mesh, local_node_map, interFace);
    output_element_fields(output_region);
    output_nodal_nodeset_fields(output_region, part_mesh, interFace);
    if (!interFace.omit_nodesets()) {
      output_nodeset_fields(output_region, merged > 0);
    }
    if (!interFace.omit_sidesets()) {
      output_sideset_fields(output_region);
    }

    for (size_t p = 0; p < part_mesh.size(); p++) {
      if (steps[p] != 0) {
        part_mesh[p]->end_state(steps[p]);
      }
    }
  }

  bool define_global_fields(Ioss::Region &output_region, const RegionVector &part_mesh,
                            const StringIdVector &variable_list)
  {
    bool error = false;
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return error;
    }
    for (const auto &pm : part_mesh) {
      Ioss::NameList fields = pm->field_describe(Ioss::Field::REDUCTION);
      for (const auto &field_name : fields) {
        if (valid_variable(field_name, 0, variable_list)) {
          Ioss::Field field = pm->get_field(field_name);
          output_region.field_add(std::move(field));
        }
      }
    }
    // Now that we have defined all fields, check `variable_list` and make
    // sure that all fields that have been explicitly specified now exist
    // on `output_region`...
    if (!variable_list.empty() && variable_list[0].first != "all") {
      // The user has specified at least one variable...
      Ioss::NameList fields = output_region.field_describe(Ioss::Field::REDUCTION);
      error                 = check_variable_mismatch("Global", variable_list, fields);
    }
    return error;
  }

  bool define_nodal_fields(const Ioss::Region &output_region, const RegionVector &part_mesh,
                           const StringIdVector &variable_list, const SystemInterface &interFace)
  {
    bool error = false;
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return error;
    }
    Ioss::NodeBlock *onb = output_region.get_node_blocks()[0];
    SMART_ASSERT(onb != nullptr);
    size_t node_count = onb->entity_count();
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      if (!interFace.convert_nodes_to_nodesets(p + 1)) {
        Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
        SMART_ASSERT(nb != nullptr);
        Ioss::NameList fields = nb->field_describe(Ioss::Field::TRANSIENT);
        for (const auto &field_name : fields) {
          if (valid_variable(field_name, 0, variable_list)) {
            Ioss::Field field = nb->get_field(field_name);
            field.reset_count(node_count);
            onb->field_add(std::move(field));
          }
        }
      }
    }
    // Now that we have defined all fields, check `variable_list` and make
    // sure that all fields that have been explicitly specified now exist
    // on `output_region`...
    if (!variable_list.empty() && variable_list[0].first != "all") {
      // The user has specified at least one variable...
      Ioss::NameList fields = onb->field_describe(Ioss::Field::REDUCTION);
      error                 = check_variable_mismatch("Nodal", variable_list, fields);
    }
    return error;
  }

  template <typename ENTITY>
  bool define_entity_fields(const StringIdVector        &variable_list,
                            const std::vector<ENTITY *> &entities, Ioss::NameList &defined_fields,
                            const std::string &type, bool check_list)
  {
    bool error             = false;
    bool subsetting_fields = !variable_list.empty() && variable_list[0].first != "all";

    for (const auto &entity : entities) {
      const auto &itr = output_input_map.find(entity);
      SMART_ASSERT(itr != output_input_map.end());
      const auto &[key, entity_inputs] = *itr;
      if (!entity_inputs.empty()) {
        int64_t count = entity->entity_count();
        for (const auto &[in_entity, offset] : entity_inputs) {
          if (in_entity != nullptr) {
            size_t         id     = in_entity->get_property("id").get_int();
            Ioss::NameList fields = in_entity->field_describe(Ioss::Field::TRANSIENT);
            for (const auto &field_name : fields) {
              if (valid_variable(field_name, id, variable_list)) {
                Ioss::Field field = in_entity->get_field(field_name);
                field.reset_count(count);
                if (!entity->field_exists(field_name)) {
                  entity->field_add(std::move(field));
                  if (subsetting_fields) {
                    defined_fields.push_back(field_name);
                  }
                }
              }
            }
          }
        }
      }
    }

    // Now that we have defined all fields, check `variable_list` and make
    // sure that all fields that have been explicitly specified now exist
    // on `output_region`...
    if (subsetting_fields && check_list) {
      // The user has specified at least one variable...
      Ioss::Utils::uniquify(defined_fields);
      error = check_variable_mismatch(type, variable_list, defined_fields);
    }
    return error;
  }

  bool define_element_fields(const Ioss::Region &output_region, const StringIdVector &variable_list)
  {
    bool error = false;
    // Element Block Fields...
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return error;
    }
    const auto &output_blocks = output_region.get_element_blocks();
    if (output_blocks.empty()) {
      return error;
    }
    Ioss::NameList defined_fields;
    return define_entity_fields(variable_list, output_blocks, defined_fields, "Element", true);
  }

  bool define_nodeset_fields(const Ioss::Region &output_region, const StringIdVector &variable_list)
  {
    bool error = false;
    // Nodeset fields...
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return error;
    }

    const auto &output_nodesets = output_region.get_nodesets();
    if (output_nodesets.empty()) {
      return error;
    }

    Ioss::NameList defined_fields;
    return define_entity_fields(variable_list, output_nodesets, defined_fields, "Nodeset", true);
  }

  void define_sideblock_fields(Ioss::SideSet *oss, const StringIdVector &variable_list,
                               Ioss::NameList &defined_fields)
  {
    // Get output sideblocks in the output sideset `oss`
    const Ioss::SideBlockContainer &osbs = oss->get_side_blocks();

    define_entity_fields(variable_list, osbs, defined_fields, "Sideblock", false);
  }

  bool define_sideset_fields(const Ioss::Region &output_region, const StringIdVector &variable_list)
  {
    bool error = false;
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return error;
    }

    bool           subsetting_fields = !variable_list.empty() && variable_list[0].first != "all";
    Ioss::NameList defined_fields;

    const auto &output_sidesets = output_region.get_sidesets();
    for (const auto &oss : output_sidesets) {
      if (output_input_map.find(oss) != output_input_map.end()) {
        define_sideblock_fields(oss, variable_list, defined_fields);
      }
    }

    // Now that we have defined all fields, check `variable_list` and make
    // sure that all fields that have been explicitly specified now exist
    // on `output_region`...
    if (subsetting_fields) {
      // The user has specified at least one variable...
      Ioss::Utils::uniquify(defined_fields);
      error = check_variable_mismatch("Sideset", variable_list, defined_fields);
    }
    return error;
  }

  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role, const std::string &prefix)
  {
    // Check for transient fields...
    auto fields = ige->field_describe(role);

    // Iterate through results fields and transfer to output
    // database...  If a prefix is specified, only transfer fields
    // whose names begin with the prefix
    for (const auto &field_name : fields) {
      if (field_name != "ids" && !oge->field_exists(field_name) &&
          Ioss::Utils::substr_equal(prefix, field_name)) {
        // If the field does not already exist, add it to the output node block
        Ioss::Field field = ige->get_field(field_name);
        oge->field_add(std::move(field));
      }
    }
  }

  bool valid_variable(const std::string &variable, size_t id, const StringIdVector &variable_list)
  {
    if (variable_list.empty() || variable_list[0].first == "all") {
      return true;
    }
    if (variable_list[0].first == "none") {
      return false;
    }

    for (const auto &[var_name, var_id] : variable_list) {
      if (var_name == variable) {
        if (id == 0 || id == var_id || var_id == 0) {
          return true;
        }
      }
    }
    return false;
  }

  bool check_variable_mismatch(const std::string &type, const StringIdVector &variable_list,
                               const Ioss::NameList &fields)
  {
    // Check all variables in `variable_list` and see if they are found in `fields`
    if (variable_list.empty() || variable_list[0].first == "all") {
      return false; // No error
    }

    bool error = false;
    for (const auto &var : variable_list) {
      if (std::find(fields.begin(), fields.end(), var.first) == std::end(fields)) {
        fmt::print(stderr, "ERROR: {} Variable '{}' was not found in the list of valid fields.\n",
                   type, var.first);
        error = true;
      }
    }
    return error;
  }

  void process_nodeset_omissions(const RegionVector &part_mesh, const Omissions &omit)
  {
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      if (!omit[p].empty()) {
        // Get the nodesets for this part and set the "omitted" property on the nodeset
        if (omit[p][0] == "ALL") {
          const Ioss::NodeSetContainer &nodesets = part_mesh[p]->get_nodesets();
          for (const auto &ns : nodesets) {
            ns->property_add(Ioss::Property(std::string("omitted"), 1));
          }
        }
        else {
          for (const auto &omitted : omit[p]) {
            Ioss::NodeSet *ns = part_mesh[p]->get_nodeset(omitted);
            if (ns != nullptr) {
              ns->property_add(Ioss::Property(std::string("omitted"), 1));
            }
          }
        }
      }
    }
  }

  void process_sideset_omissions(const RegionVector &part_mesh, const Omissions &omit)
  {
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      if (!omit[p].empty()) {
        // Get the sidesets for this part and set the "omitted" property on the sideset
        if (omit[p][0] == "ALL") {
          const Ioss::SideSetContainer &sidesets = part_mesh[p]->get_sidesets();
          for (const auto &ss : sidesets) {
            ss->property_add(Ioss::Property(std::string("omitted"), 1));
          }
        }
        else {
          for (const auto &omitted : omit[p]) {
            Ioss::SideSet *ss = part_mesh[p]->get_sideset(omitted);
            if (ss != nullptr) {
              ss->property_add(Ioss::Property(std::string("omitted"), 1));
            }
          }
        }
      }
    }
  }

  void process_assembly_omissions(const RegionVector &part_mesh, const Omissions &omit)
  {
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      if (!omit[p].empty()) {
        // Get the assemblies for this part and set the "omitted" property on the assembly
        if (omit[p][0] == "ALL") {
          const auto &assemblies = part_mesh[p]->get_assemblies();
          for (const auto &as : assemblies) {
            as->property_add(Ioss::Property(std::string("omitted"), 1));
          }
        }
        else {
          for (const auto &omitted : omit[p]) {
            auto *as = part_mesh[p]->get_assembly(omitted);
            if (as != nullptr) {
              as->property_add(Ioss::Property(std::string("omitted"), 1));
            }
          }
        }
      }
    }
  }

  template <typename INT>
  void check_for_duplicate_nodeset_nodes(Ioss::Region           &output_region,
                                         const std::vector<INT> &local_node_map)
  {
    const auto &output_nodesets = output_region.get_nodesets();
    for (const auto &ons : output_nodesets) {
      const auto &itr = output_input_map.find(ons);
      SMART_ASSERT(itr != output_input_map.end());
      const auto &[key, ons_inputs] = *itr;
      if (ons_inputs.size() >= 2) {
        int64_t count = 0;
        for (const auto &[ins, offset] : ons_inputs) {
          if (ins != nullptr) {
            count += ins->entity_count();
          }
        }

        std::vector<INT> nodelist(count);
        int              found = 0;
        for (const auto &[ins, offset] : ons_inputs) {
          if (ins != nullptr) {
            ++found;
            ins->get_field_data("ids", &nodelist[offset], -1);

            auto *input_region = dynamic_cast<const Ioss::Region *>(ins->contained_in());
            Ioss::Utils::check_dynamic_cast(input_region);
            size_t node_offset = input_region->get_property("node_offset").get_int();
            for (int64_t i = 0; i < ins->entity_count(); i++) {
              size_t loc_node = input_region->node_global_to_local(nodelist[offset + i], true) - 1;
              auto   gpos     = local_node_map[node_offset + loc_node];
              if (gpos >= 0) {
                nodelist[offset + i] = gpos + 1;
              }
            }
          }
        }
        if (found >= 2) {
          auto size_pre  = nodelist.size();
          auto size_post = size_pre;
          {
            auto tmp_nodelist = nodelist;
            Ioss::Utils::uniquify(tmp_nodelist);
            size_post = tmp_nodelist.size();
          }
          if (size_pre != size_post) {
            reset_entity_count(ons, size_post);

            // Create a map for `ons` that maps the combined ins positions into the
            // output position accounting for eliminating duplicate nodes.
            std::vector<std::pair<INT, INT>> ids_pos;
            ids_pos.reserve(size_pre);
            for (auto [i, id] : Ioss::enumerate(nodelist)) {
              ids_pos.emplace_back(id, i);
            }
            std::sort(ids_pos.begin(), ids_pos.end());

            auto new_size = unique(ids_pos);
            SMART_ASSERT(new_size == (size_t)ons->entity_count())(new_size)(ons->entity_count());
            SMART_ASSERT(new_size == size_post)(new_size)(size_post);

            auto &map_vector = nodeset_in_out_map[ons];
            SMART_ASSERT(map_vector.empty())(map_vector.size());
            map_vector.reserve(new_size);
            for (const auto &[id, pos] : ids_pos) {
              map_vector.push_back(pos);
            }
            SMART_ASSERT(map_vector.size() == (size_t)ons->entity_count())
            (map_vector.size())(ons->entity_count());
          }
        }
      }
    }
  }
} // namespace
