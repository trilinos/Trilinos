// Copyright(C) 1999-2024 National Technology & Engineering Solutions
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

#include "add_to_log.h"
#include "fmt/ostream.h"
#include "time_stamp.h"

#include <exodusII.h>

#include <Ionit_Initializer.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Transform.h>

#include "EJ_CodeTypes.h"
#include "EJ_SystemInterface.h"
#include "EJ_Version.h"
#include "EJ_mapping.h"
#include "EJ_match_xyz.h"
#include "EJ_vector3d.h"

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif

namespace {
  bool valid_variable(const std::string &variable, size_t id, const StringIdVector &variable_list);
  void define_global_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                            const StringIdVector &variable_list);
  void define_nodal_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                           const StringIdVector &variable_list, SystemInterface &interFace);
  void define_element_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                             const StringIdVector &variable_list);
  void define_nset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                          const StringIdVector &variable_list);
  void define_sset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                          const StringIdVector &variable_list);
  void define_nodal_nodeset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                                   const StringIdVector &variable_list, SystemInterface &interFace);

  template <typename INT>
  void output_nodeblock(Ioss::Region &output_region, RegionVector &part_mesh,
                        const std::vector<INT> &local_node_map, std::vector<INT> &global_node_map);
  template <typename INT>
  void output_elementblock(Ioss::Region &output_region, RegionVector &part_mesh,
                           const std::vector<INT> &local_node_map,
                           const std::vector<INT> &local_element_map, bool ignore_element_ids);
  template <typename INT>
  void output_nodeset(Ioss::Region &output_region, RegionVector &part_mesh,
                      const std::vector<INT> &local_node_map);
  template <typename INT>
  void output_sideset(Ioss::Region &output_region, RegionVector &part_mesh,
                      const std::vector<INT> &local_element_map);
  template <typename INT>
  void output_nodal_nodeset(Ioss::Region &output_region, RegionVector &part_mesh,
                            SystemInterface &interFace, const std::vector<INT> &local_node_map);
  template <typename INT>
  void output_transient_state(Ioss::Region &output_region, RegionVector &part_mesh, double time,
                              const std::vector<INT> &local_node_map, SystemInterface &interFace);
  void process_nset_omissions(RegionVector &part_mesh, const Omissions &omit);
  void process_sset_omissions(RegionVector &part_mesh, const Omissions &omit);
  void process_assembly_omissions(RegionVector &part_mesh, const Omissions &omit);

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

} // namespace

namespace {
  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region,
                             bool create_assemblies, bool debug);
  void transfer_assembly(const Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void create_nodal_nodeset(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_fields(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role, const std::string &prefix = "");

  void transfer_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                           Ioss::Field::RoleType role, const std::string &prefix = "",
                           bool transfer_connectivity = true);

  void transfer_field_data_internal(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                                    const std::string &field_name);
} // namespace

std::string tsFormat = "[%H:%M:%S] ";

// prototypes

template <typename INT>
double ejoin(SystemInterface &interFace, std::vector<Ioss::Region *> &part_mesh, INT dummy);

unsigned int debug_level = 0;

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
    std::vector<Ioss::Region *>     part_mesh(interFace.inputFiles_.size());
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

      vector3d offset = interFace.offset();
      if (p > 0 && (offset.x != 0.0 || offset.y != 0.0 || offset.z != 0.0)) {
        Ioss::NodeBlock *nb        = part_mesh[p]->get_node_blocks()[0];
        Ioss::Field      coord     = nb->get_field("mesh_model_coordinates");
        auto            *transform = Ioss::Transform::create("offset3D");
        assert(transform != nullptr);
        std::vector<double> values(3);
        values[0] = offset.x * p;
        values[1] = offset.y * p;
        values[2] = offset.z * p;
        transform->set_properties("offset", values);
        coord.add_transform(transform);
        nb->field_erase("mesh_model_coordinates");
        nb->field_add(coord);
      }
    }

    process_nset_omissions(part_mesh, interFace.nset_omissions());
    process_sset_omissions(part_mesh, interFace.sset_omissions());
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

template <typename INT>
double ejoin(SystemInterface &interFace, std::vector<Ioss::Region *> &part_mesh, INT /*dummy*/)
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

  if (interFace.compression_level() > 0 || interFace.szip()) {
    properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    properties.add(Ioss::Property("COMPRESSION_LEVEL", interFace.compression_level()));
    properties.add(Ioss::Property("COMPRESSION_SHUFFLE", true));
    if (interFace.szip()) {
      properties.add(Ioss::Property("COMPRESSION_METHOD", "szip"));
    }
    else if (interFace.zlib()) {
      properties.add(Ioss::Property("COMPRESSION_METHOD", "zlib"));
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

  INT node_offset    = 0;
  INT element_offset = 0;
  for (auto &pm : part_mesh) {
    pm->property_add(Ioss::Property("node_offset", node_offset));
    pm->property_add(Ioss::Property("element_offset", element_offset));
    INT local_node_count = pm->get_property("node_count").get_int();
    INT local_elem_count = pm->get_property("element_count").get_int();
    node_offset += local_node_count;
    element_offset += local_elem_count;
  }

  INT              node_count = node_offset; // Sum of nodes in part meshes.
  std::vector<INT> local_node_map(node_count);
  std::vector<INT> global_node_map;

  // This is the map from local element position to global element
  // position (0-based). If there are no element block omissions, then
  // the map is simply [0..number_elements). If there are omissions,
  // Then local_element_map[j] will be -1 for an omitted element.
  // If not omitted, local_element_map[part_offset+j] gives the global
  // position of local element j in the current part.
  std::vector<INT> local_element_map(element_offset);
  build_local_element_map(part_mesh, local_element_map);

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

  // Add element blocks, nodesets, sidesets
  for (size_t p = 0; p < part_count; p++) {
    transfer_elementblock(*part_mesh[p], output_region, interFace.create_assemblies(), false);
    if (interFace.convert_nodes_to_nodesets(p + 1)) {
      create_nodal_nodeset(*part_mesh[p], output_region, false);
    }
    if (!interFace.omit_nodesets()) {
      transfer_nodesets(*part_mesh[p], output_region, false);
    }
    if (!interFace.omit_sidesets()) {
      transfer_sidesets(*part_mesh[p], output_region, false);
    }
    if (!interFace.omit_assemblies()) {
      transfer_assembly(*part_mesh[p], output_region, false);
    }
  }

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
                      interFace.ignore_element_ids());
  output_nodal_nodeset(output_region, part_mesh, interFace, local_node_map);

  if (!interFace.omit_nodesets()) {
    output_nodeset(output_region, part_mesh, local_node_map);
  }
  if (!interFace.omit_sidesets()) {
    output_sideset(output_region, part_mesh, local_element_map);
  }

  output_region.end_mode(Ioss::STATE_MODEL);

  // ####################TRANSIENT DATA SECTION###########################
  // ***********************************************************************
  // 9. Get Variable Information and names

  if (debug_level & 1) {
    fmt::print(stderr, "{}", time_stamp(tsFormat));
  }

  output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

  define_global_fields(output_region, part_mesh, interFace.global_var_names());

  define_nodal_fields(output_region, part_mesh, interFace.node_var_names(), interFace);
  define_nodal_nodeset_fields(output_region, part_mesh, interFace.node_var_names(), interFace);

  define_element_fields(output_region, part_mesh, interFace.elem_var_names());

  if (!interFace.omit_nodesets()) {
    define_nset_fields(output_region, part_mesh, interFace.nset_var_names());
  }
  if (!interFace.omit_sidesets()) {
    define_sset_fields(output_region, part_mesh, interFace.sset_var_names());
  }

  output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

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
    output_transient_state(output_region, part_mesh, global_times[step], local_node_map, interFace);
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

  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region,
                             bool create_assemblies, bool debug)
  {
    static int         used_blocks = 0;
    const std::string &prefix      = region.name();

    Ioss::Assembly *assem = nullptr;

    const Ioss::ElementBlockContainer &ebs = region.get_element_blocks();
    for (const auto &eb : ebs) {
      if (!entity_is_omitted(eb)) {
        std::string name = eb->name();
        if (output_region.get_element_block(name) != nullptr) {
          name = prefix + "_" + eb->name();
          if (output_region.get_element_block(name) != nullptr) {
            fmt::print(stderr, "ERROR: Duplicate element blocks named '{}'\n", name);
            exit(EXIT_FAILURE);
          }
        }
        eb->property_add(Ioss::Property("name_in_output", name));

        if (debug) {
          fmt::print(stderr, "{}, ", name);
        }
        std::string type     = eb->topology()->name();
        size_t      num_elem = eb->entity_count();

        if (num_elem > 0) {
          auto *ebn = new Ioss::ElementBlock(output_region.get_database(), name, type, num_elem);
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
            assert(entity != nullptr);
            asn->add(entity);
          }
        }
      }
    }
  }

  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const std::string &prefix = region.name();

    const Ioss::SideSetContainer &fss = region.get_sidesets();
    for (const auto &fs : fss) {
      if (!entity_is_omitted(fs)) {
        std::string name = fs->name();
        if (output_region.get_sideset(name) != nullptr) {
          name = prefix + "_" + fs->name();
          if (output_region.get_sideset(name) != nullptr) {
            fmt::print(stderr, "ERROR: Duplicate side sets named '{}'\n", name);
            exit(EXIT_FAILURE);
          }
        }
        fs->property_add(Ioss::Property("name_in_output", name));
        if (debug) {
          fmt::print(stderr, "{}, ", name);
        }
        auto *surf = new Ioss::SideSet(output_region.get_database(), name);
        set_id(fs, surf);

        const Ioss::SideBlockContainer &fbs = fs->get_side_blocks();
        for (const auto &fb : fbs) {
          const std::string &fbname = prefix + "_" + fb->name();
          if (debug) {
            fmt::print(stderr, "{}, ", fbname);
          }
          fb->property_add(Ioss::Property("name_in_output", fbname));
          std::string fbtype   = fb->topology()->name();
          std::string partype  = fb->parent_element_topology()->name();
          size_t      num_side = fb->entity_count();

          auto *block =
              new Ioss::SideBlock(output_region.get_database(), fbname, fbtype, partype, num_side);
          surf->add(block);
        }
        output_region.add(surf);
      }
    }
  }

  // Create a nodeset on the output region consisting of all the nodes
  // in the input region.
  void create_nodal_nodeset(Ioss::Region &region, Ioss::Region &output_region, bool debug)
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
  }

  // Output the bulk data for a nodeset on the output region
  // consisting of all the nodes in the input region.
  template <typename INT>
  void output_nodal_nodeset(Ioss::Region &output_region, RegionVector &part_mesh,
                            SystemInterface &interFace, const std::vector<INT> &local_node_map)
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

  void define_nodal_nodeset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                                   const StringIdVector &variable_list, SystemInterface &interFace)
  {
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

  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    const std::string &prefix = region.name();

    const Ioss::NodeSetContainer &nss = region.get_nodesets();
    for (const auto &ns : nss) {
      if (!entity_is_omitted(ns)) {
        std::string name = ns->name();
        if (output_region.get_nodeset(name) != nullptr) {
          name = prefix + "_" + ns->name();
          if (output_region.get_nodeset(name) != nullptr) {
            fmt::print(stderr, "ERROR: Duplicate node sets named '{}'\n", name);
            exit(EXIT_FAILURE);
          }
        }
        ns->property_add(Ioss::Property("name_in_output", name));
        if (debug) {
          fmt::print(stderr, "{}, ", name);
        }
        size_t count    = ns->entity_count();
        auto  *node_set = new Ioss::NodeSet(output_region.get_database(), name, count);
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
  void output_nodeblock(Ioss::Region &output_region, RegionVector &part_mesh,
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
  void output_elementblock(Ioss::Region &output_region, RegionVector &part_mesh,
                           const std::vector<INT> &local_node_map,
                           const std::vector<INT> &local_element_map, bool ignore_element_ids)
  {

    const Ioss::ElementBlockContainer &ebs = output_region.get_element_blocks();

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
    for (const auto &pm : part_mesh) {
      const Ioss::ElementBlockContainer &iebs        = pm->get_element_blocks();
      size_t                             node_offset = pm->get_property("node_offset").get_int();

      for (const auto &ieb : iebs) {
        if (entity_is_omitted(ieb)) {
          continue;
        }
        std::string         name = pm->name() + "_" + ieb->name();
        Ioss::ElementBlock *oeb  = output_region.get_element_block(name);
        if (oeb == nullptr) {
          name = ieb->name();
          oeb  = output_region.get_element_block(name);
        }
        if (oeb != nullptr) {
          std::vector<INT> connectivity;
          ieb->get_field_data("connectivity_raw", connectivity);

          SMART_ASSERT(ieb->entity_count() == oeb->entity_count());
          for (auto &node : connectivity) {
            // connectivity is in part-local node ids [1..num_node]
            // loc_node = the position of node in the local [0..num_node)
            // local_node_map[node_offset+loc_node] gives the position of this node in the global
            // list
            size_t loc_node = node - 1;
            SMART_ASSERT(node_offset + loc_node < local_node_map.size());
            auto gpos = local_node_map[node_offset + loc_node];
            if (gpos >= 0) {
              node = gpos + 1;
            }
          }
          oeb->put_field_data("connectivity_raw", connectivity);
          transfer_field_data(ieb, oeb, Ioss::Field::ATTRIBUTE);
        }
      }
    }
  }

  template <typename INT>
  void output_nodeset(Ioss::Region &output_region, RegionVector &part_mesh,
                      const std::vector<INT> &local_node_map)
  {
    if (output_region.get_nodesets().empty()) {
      return;
    }

    for (const auto &pm : part_mesh) {
      size_t                        node_offset = pm->get_property("node_offset").get_int();
      const Ioss::NodeSetContainer &ins         = pm->get_nodesets();
      for (const auto &in : ins) {
        if (!entity_is_omitted(in)) {
          std::vector<INT> nodelist;
          in->get_field_data("ids", nodelist);

          std::string    name = pm->name() + "_" + in->name();
          Ioss::NodeSet *ons  = output_region.get_nodeset(name);
          if (ons == nullptr) {
            name = in->name();
            ons  = output_region.get_nodeset(name);
          }
          SMART_ASSERT(ons != nullptr)(name);
          SMART_ASSERT(in->entity_count() == ons->entity_count());

          // This needs to make sure that the nodelist comes back as local id (1..numnodes)
          for (auto &node : nodelist) {
            size_t loc_node = pm->node_global_to_local(node, true) - 1;
            auto   gpos     = local_node_map[node_offset + loc_node];
            if (gpos >= 0) {
              node = gpos + 1;
            }
          }
          ons->put_field_data("ids_raw", nodelist);

          std::vector<double> df;
          in->get_field_data("distribution_factors", df);
          ons->put_field_data("distribution_factors", df);
        }
      }
    }
  }

  template <typename INT>
  void output_sideset(Ioss::Region &output_region, RegionVector &part_mesh,
                      const std::vector<INT> &local_element_map)
  {
    const Ioss::SideSetContainer &os = output_region.get_sidesets();

    Ioss::SideBlockContainer out_eb;
    // Put all output side blocks in the same list...
    for (const auto &oss : os) {
      const Ioss::SideBlockContainer &obs = oss->get_side_blocks();
      std::copy(obs.begin(), obs.end(), std::back_inserter(out_eb));
    }

    // Assuming (with checks) that the output side blocks will be
    // iterated in same order as input side blocks...
    auto II = out_eb.begin();

    for (const auto &pm : part_mesh) {
      size_t element_offset = pm->get_property("element_offset").get_int();

      const Ioss::SideSetContainer &is = pm->get_sidesets();
      for (const auto &iss : is) {
        if (!entity_is_omitted(iss)) {
          const Ioss::SideBlockContainer &ebs = iss->get_side_blocks();

          for (const auto &eb : ebs) {
            SMART_ASSERT((eb->name() == (*II)->name()) ||
                         (pm->name() + "_" + eb->name() == (*II)->name()))
            (eb->name())((*II)->name());
            SMART_ASSERT(eb->entity_count() == (*II)->entity_count());
            std::vector<INT> elem_side_list;
            eb->get_field_data("element_side_raw", elem_side_list);

            // The 'elem_side_list' contains
            // (local_element_position,side_ordinal) pairs. The
            // 'local_element_position' is 1-based offset in the
            // current part.  Need to map to its location in the
            // output region...
            for (size_t i = 0; i < elem_side_list.size();
                 i += 2) { // just get the elem part of the pair...
              size_t local_position = elem_side_list[i] - 1;
              auto   gpos           = local_element_map[element_offset + local_position];
              SMART_ASSERT(gpos >= 0)(gpos)(i); // Inactive elements should be filtered by Ioss
              elem_side_list[i] = gpos + 1;
            }
            (*II)->put_field_data("element_side_raw", elem_side_list);
            ++II;
          }
        }
      }
    }
  }

  void output_globals(Ioss::Region &output_region, RegionVector &part_mesh)
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
  void output_nodal(Ioss::Region &output_region, RegionVector &part_mesh,
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

  void output_nodal_nodeset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
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

  void output_element(Ioss::Region &output_region, RegionVector &part_mesh)
  {
    for (const auto &pm : part_mesh) {
      const Ioss::ElementBlockContainer &iebs = pm->get_element_blocks();
      for (const auto &ieb : iebs) {
        if (!entity_is_omitted(ieb)) {
          std::string         name = pm->name() + "_" + ieb->name();
          Ioss::ElementBlock *oeb  = output_region.get_element_block(name);
          if (oeb == nullptr) {
            name = ieb->name();
            oeb  = output_region.get_element_block(name);
          }
          if (oeb != nullptr) {
            Ioss::NameList fields = ieb->field_describe(Ioss::Field::TRANSIENT);
            for (const auto &field : fields) {
              if (oeb->field_exists(field)) {
                transfer_field_data_internal(ieb, oeb, field);
              }
            }
          }
        }
      }
    }
  }

  void output_nset(Ioss::Region &output_region, RegionVector &part_mesh)
  {
    if (output_region.get_nodesets().empty()) {
      return;
    }

    for (const auto &pm : part_mesh) {
      const Ioss::NodeSetContainer &ins = pm->get_nodesets();
      for (const auto &in : ins) {
        if (!entity_is_omitted(in)) {
          std::string    name = pm->name() + "_" + in->name();
          Ioss::NodeSet *ons  = output_region.get_nodeset(name);
          if (ons == nullptr) {
            name = in->name();
            ons  = output_region.get_nodeset(name);
          }
          SMART_ASSERT(ons != nullptr)(name);

          Ioss::NameList fields = in->field_describe(Ioss::Field::TRANSIENT);
          for (const auto &field : fields) {
            if (ons->field_exists(field)) {
              transfer_field_data_internal(in, ons, field);
            }
          }
        }
      }
    }
  }

  void output_sset(Ioss::Region &output_region, RegionVector &part_mesh)
  {
    const Ioss::SideSetContainer &os = output_region.get_sidesets();
    if (os.empty()) {
      return;
    }

    Ioss::SideBlockContainer out_eb;
    // Put all output side blocks in the same list...
    for (const auto &oss : os) {
      const Ioss::SideBlockContainer &obs = oss->get_side_blocks();
      std::copy(obs.begin(), obs.end(), std::back_inserter(out_eb));
    }

    // Assuming (with checks) that the output side blocks will be
    // iterated in same order as input side blocks...
    auto II = out_eb.begin();

    for (const auto &pm : part_mesh) {
      const Ioss::SideSetContainer &is = pm->get_sidesets();
      for (const auto &iss : is) {
        if (!entity_is_omitted(iss)) {
          const Ioss::SideBlockContainer &ebs = iss->get_side_blocks();
          for (const auto &eb : ebs) {
            SMART_ASSERT((pm->name() + "_" + eb->name() == (*II)->name()) ||
                         (eb->name() == (*II)->name()));
            Ioss::NameList fields = eb->field_describe(Ioss::Field::TRANSIENT);
            for (const auto &field : fields) {
              if ((*II)->field_exists(field)) {
                transfer_field_data_internal(eb, *II, field);
              }
            }
            ++II;
          }
        }
      }
    }
  }

  template <typename INT>
  void output_transient_state(Ioss::Region &output_region, RegionVector &part_mesh, double time,
                              const std::vector<INT> &local_node_map, SystemInterface &interFace)
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

    output_globals(output_region, part_mesh);
    output_nodal(output_region, part_mesh, local_node_map, interFace);
    output_element(output_region, part_mesh);
    output_nodal_nodeset_fields(output_region, part_mesh, interFace);
    if (!interFace.omit_nodesets()) {
      output_nset(output_region, part_mesh);
    }
    if (!interFace.omit_sidesets()) {
      output_sset(output_region, part_mesh);
    }

    for (size_t p = 0; p < part_mesh.size(); p++) {
      if (steps[p] != 0) {
        part_mesh[p]->end_state(steps[p]);
      }
    }
  }

  void transfer_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                           Ioss::Field::RoleType role, const std::string &prefix,
                           bool transfer_connectivity)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields = ige->field_describe(role);

    // Complication here is that if the 'role' is 'Ioss::Field::MESH',
    // then the 'ids' field must be transferred first...
    if (role == Ioss::Field::MESH) {
      for (const auto &field_name : state_fields) {
        if (oge->field_exists(field_name)) {
          if (field_name == "ids") {
            transfer_field_data_internal(ige, oge, field_name);
            break;
          }
        }
      }
    }

    for (const auto &field_name : state_fields) {
      // All of the 'Ioss::EntityBlock' derived classes have a
      // 'connectivity' field, but it is only interesting on the
      // Ioss::ElementBlock class. On the other classes, it just
      // generates overhead...
      if (!transfer_connectivity && field_name == "connectivity") {
        continue;
      }

      if (field_name != "ids" && Ioss::Utils::substr_equal(prefix, field_name)) {
        if (oge->field_exists(field_name)) {
          transfer_field_data_internal(ige, oge, field_name);
        }
      }
    }
  }

  void transfer_field_data_internal(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                                    const std::string &field_name)
  {
    static std::vector<double> data;

    assert(ige->get_field(field_name).get_size() == oge->get_field(field_name).get_size());

    ige->get_field_data(field_name, data);
    oge->put_field_data(field_name, data);
  }

  void define_global_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                            const StringIdVector &variable_list)
  {
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return;
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
  }

  void define_nodal_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                           const StringIdVector &variable_list, SystemInterface &interFace)
  {
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return;
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
  }

  void define_element_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                             const StringIdVector &variable_list)
  {
    // Element Block Fields...
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return;
    }
    for (const auto &pm : part_mesh) {
      const Ioss::ElementBlockContainer &iebs = pm->get_element_blocks();
      for (const auto &ieb : iebs) {
        if (!entity_is_omitted(ieb)) {
          std::string         name = pm->name() + "_" + ieb->name();
          Ioss::ElementBlock *oeb  = output_region.get_element_block(name);
          if (oeb == nullptr) {
            name = ieb->name();
            oeb  = output_region.get_element_block(name);
          }
          if (oeb != nullptr) {
            size_t         id     = oeb->get_property("id").get_int();
            Ioss::NameList fields = ieb->field_describe(Ioss::Field::TRANSIENT);
            for (const auto &field_name : fields) {
              if (valid_variable(field_name, id, variable_list)) {
                Ioss::Field field = ieb->get_field(field_name);
                oeb->field_add(std::move(field));
              }
            }
          }
        }
      }
    }
  }

  void define_nset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                          const StringIdVector &variable_list)
  {
    // Nodeset fields...
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return;
    }

    for (const auto &pm : part_mesh) {
      const Ioss::NodeSetContainer &ins = pm->get_nodesets();
      for (const auto &in : ins) {
        if (!entity_is_omitted(in)) {
          std::string    name = pm->name() + "_" + in->name();
          Ioss::NodeSet *ons  = output_region.get_nodeset(name);
          if (ons == nullptr) {
            name = in->name();
            ons  = output_region.get_nodeset(name);
          }
          SMART_ASSERT(ons != nullptr)(name);

          size_t         id     = in->get_property("id").get_int();
          Ioss::NameList fields = in->field_describe(Ioss::Field::TRANSIENT);
          for (const auto &field_name : fields) {
            if (valid_variable(field_name, id, variable_list)) {
              Ioss::Field field = in->get_field(field_name);
              ons->field_add(std::move(field));
            }
          }
        }
      }
    }
  }

  void define_sset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
                          const StringIdVector &variable_list)
  {
    if (!variable_list.empty() && variable_list[0].first == "none") {
      return;
    }
    const auto &os = output_region.get_sidesets();

    Ioss::SideBlockContainer out_eb;
    // Put all output side blocks in the same list...
    for (const auto &oss : os) {
      const auto &obs = oss->get_side_blocks();
      std::copy(obs.begin(), obs.end(), std::back_inserter(out_eb));
    }

    // Assuming (with checks) that the output side blocks will be
    // iterated in same order as input side blocks...
    auto II = out_eb.begin();

    for (const auto &pm : part_mesh) {
      const auto &is = pm->get_sidesets();
      for (const auto &iss : is) {
        if (!entity_is_omitted(iss)) {
          size_t      id  = iss->get_property("id").get_int();
          const auto &ebs = iss->get_side_blocks();
          for (const auto &eb : ebs) {
            SMART_ASSERT((pm->name() + "_" + eb->name() == (*II)->name()) ||
                         (eb->name() == (*II)->name()));
            auto fields = eb->field_describe(Ioss::Field::TRANSIENT);
            for (const auto &field_name : fields) {
              if (valid_variable(field_name, id, variable_list)) {
                Ioss::Field field = eb->get_field(field_name);
                (*II)->field_add(std::move(field));
              }
            }
            ++II;
          }
        }
      }
    }
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

    for (const auto &var : variable_list) {
      if (var.first == variable) {
        if (id == 0 || id == var.second || var.second == 0) {
          return true;
        }
      }
    }
    return false;
  }

  void process_nset_omissions(RegionVector &part_mesh, const Omissions &omit)
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

  void process_sset_omissions(RegionVector &part_mesh, const Omissions &omit)
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

  void process_assembly_omissions(RegionVector &part_mesh, const Omissions &omit)
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
} // namespace
