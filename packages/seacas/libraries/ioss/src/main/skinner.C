// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>

#include "Ionit_Initializer.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DataSize.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_EntityBlock.h"
#include "Ioss_FaceGenerator.h"
#include "Ioss_Field.h"
#include "Ioss_FileInfo.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_State.h"
#include "Ioss_Utils.h"
#include "robin_hash.h"
#include "skinner_interface.h"

// ========================================================================
namespace {
  template <typename INT> void skinner(Skinner::Interface &interFace, INT /*dummy*/);
  std::string                  codename;
  std::string                  version = "1.02";

  void transfer_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                           Ioss::Field::RoleType role, const Ioss::IntVector &ref_nodes);

  void transfer_field_data(Ioss::EntityBlock *ige, Ioss::EntityBlock *oge,
                           Ioss::Field::RoleType          role,
                           const std::vector<Ioss::Face> &boundary_faces);

  void output_table(const Ioss::ElementBlockContainer              &ebs,
                    std::map<std::string, std::vector<Ioss::Face>> &boundary_faces)
  {
    // Get maximum name and face_count length...
    size_t max_name = std::string("Block Name").length();
    size_t max_face = std::string("Face Count").length();
    for (auto &eb : ebs) {
      const std::string &name = eb->name();
      if (name.length() > max_name) {
        max_name = name.length();
      }
      size_t face_width = Ioss::Utils::number_width(boundary_faces[name].size());
      max_face          = face_width > max_face ? face_width : max_face;
    }
    max_name += 4; // Padding
    max_face += 4;

    fmt::print("\t+{2:-^{0}}+{2:-^{1}}+\n", max_name, max_face, "");
    fmt::print("\t|{2:^{0}}|{3:^{1}}|\n", max_name, max_face, "Block Name", "Face Count");
    fmt::print("\t+{2:-^{0}}+{2:-^{1}}+\n", max_name, max_face, "");
    for (auto &eb : ebs) {
      const std::string &name = eb->name();
      fmt::print("\t|{2:^{0}}|{3:{1}}  |\n", max_name, max_face - 2, name,
                 fmt::group_digits(boundary_faces[name].size()));
    }
    fmt::print("\t+{2:-^{0}}+{2:-^{1}}+\n", max_name, max_face, "");
  }

  Ioss::IntVector get_selected_steps(Ioss::Region &region, const Skinner::Interface &options)
  {
    // This routine checks all steps of the input database and selects those which
    // meet the requirements specified in `options`.  The returned (1-based) vector will have a
    // value of `1` if the step is to be output and `0` if skipped.
    int             step_count = (int)region.get_property("state_count").get_int();
    Ioss::IntVector selected_steps(step_count + 1);

    // If user specified a list of times to transfer to output database,
    // process the list and find the times on the input database that are
    // closest to the times in the list.
    if (!options.selected_times.empty()) {
      int selected_step = 0;
      for (auto time : options.selected_times) {
        double diff = std::numeric_limits<double>::max();
        for (int step = 1; step <= step_count; step++) {
          double db_time  = region.get_state_time(step);
          double cur_diff = std::abs(db_time - time);
          if (cur_diff < diff) {
            diff          = std::abs(db_time - time);
            selected_step = step;
          }
        }
        if (selected_step > 0) {
          selected_steps[selected_step] = 1;
        }
      }
    }
    else {
      // User did not select specific times to be output...
      // Just select them all
      for (int i = 1; i <= step_count; i++) {
        selected_steps[i] = 1;
      }
    }

    // Now, filter by min and max time...
    for (int istep = 1; istep <= step_count; istep++) {
      double time = region.get_state_time(istep);
      if (time < options.minimum_time) {
        selected_steps[istep] = 0;
      }
      if (time > options.maximum_time) {
        selected_steps[istep] = 0;
      }
    }
    return selected_steps;
  }
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif
  Ioss::ParallelUtils pu{};
  int                 my_rank = pu.parallel_rank();

  codename = Ioss::FileInfo(argv[0]).basename();

  Skinner::Interface interFace;
  bool               success = interFace.parse_options(argc, argv);
  if (!success) {
    return EXIT_FAILURE;
  }

  Ioss::Init::Initializer io;

  if (my_rank == 0) {
    fmt::print("\nInput:    '{}', Type: {}\n", interFace.input_filename(), interFace.input_type());
    if (!interFace.no_output()) {
      fmt::print("Output:   '{}', Type: {}\n", interFace.output_filename(),
                 interFace.output_type());
    }
  }

  double begin = Ioss::Utils::timer();
  try {
    if (interFace.ints_64_bit()) {
      skinner(interFace, static_cast<int64_t>(0));
    }
    else {
      skinner(interFace, 0);
    }
  }
  catch (std::exception &e) {
    fmt::print(stderr, "\n{}\n\nskinner terminated due to exception\n", e.what());
    exit(EXIT_FAILURE);
  }
  pu.barrier();
  double end = Ioss::Utils::timer();

  if (my_rank == 0) {
    fmt::print("\n\tTotal Execution Time = {:.4} seconds\n", end - begin);
    fmt::print("\n{} execution successful.\n\n", codename);
  }
  return EXIT_SUCCESS;
}

namespace {
  Ioss::PropertyManager set_properties(Skinner::Interface &interFace);

  template <typename INT> void skinner(Skinner::Interface &interFace, INT /*dummy*/)
  {
    std::string inpfile    = interFace.input_filename();
    std::string input_type = interFace.input_type();

    Ioss::PropertyManager properties = set_properties(interFace);

    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART,
                                                    Ioss::ParallelUtils::comm_world(), properties);
    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    if (interFace.ints_64_bit()) {
      dbi->set_int_byte_size_api(Ioss::USE_INT64_API);
    }

    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");
    int          my_rank = region.get_database()->util().parallel_rank();

    if (my_rank == 0) {
      region.output_summary(std::cerr, false);
    }

    // Generate the faces...
    Ioss::FaceGenerator face_generator(region);
    face_generator.generate_faces((INT)0, true);

    // Get vector of all boundary faces which will be output as the skin...

    std::map<std::string, std::vector<Ioss::Face>> boundary_faces;
    const Ioss::ElementBlockContainer             &ebs = region.get_element_blocks();
    for (auto &eb : ebs) {
      const std::string &name     = eb->name();
      auto              &boundary = boundary_faces[name];
      auto              &faces    = face_generator.faces(name);
      for (auto &face : faces) {
        if (face.element_count() == 1) {
          boundary.push_back(face);
        }
      }
    }
    output_table(ebs, boundary_faces);

    if (interFace.no_output()) {
      return;
    }

    // Iterate the boundary faces and determine which nodes are referenced...
    size_t          my_node_count = region.get_property("node_count").get_int();
    Ioss::IntVector ref_nodes(my_node_count);
    for (const auto &boundaries : boundary_faces) {
      for (const auto &face : boundaries.second) {
        for (auto &gnode : face.connectivity_) {
          if (gnode > 0) {
            auto node       = region.get_database()->node_global_to_local(gnode, true) - 1;
            ref_nodes[node] = 1;
          }
        }
      }
    }

    size_t ref_count = std::accumulate(ref_nodes.begin(), ref_nodes.end(), 0);

    // Map ids from total mesh down to skin mesh...
    // Also get coordinates...
    std::vector<double> coord_in;
    Ioss::NodeBlock    *nb = region.get_node_blocks()[0];
    nb->get_field_data("mesh_model_coordinates", coord_in);

    Ioss::IntVector ids;
    nb->get_field_data("ids", ids);

    Ioss::IntVector owner;
    nb->get_field_data("owning_processor", owner);

    Ioss::IntVector     ref_ids(ref_count);
    std::vector<double> coord_out(3 * ref_count);
    Ioss::IntVector     owner_out(ref_count);

    size_t j = 0;
    for (size_t i = 0; i < ref_nodes.size(); i++) {
      if (ref_nodes[i] == 1) {
        coord_out[3 * j + 0] = coord_in[3 * i + 0];
        coord_out[3 * j + 1] = coord_in[3 * i + 1];
        coord_out[3 * j + 2] = coord_in[3 * i + 2];
        owner_out[j]         = owner[i];
        ref_ids[j++]         = ids[i];
      }
    }
    Ioss::Utils::clear(coord_in);
    Ioss::Utils::clear(ids);
    Ioss::Utils::clear(owner);

    std::string       file = interFace.output_filename();
    std::string       type = interFace.output_type();
    Ioss::DatabaseIO *dbo  = Ioss::IOFactory::create(type, file, Ioss::WRITE_RESTART,
                                                     Ioss::ParallelUtils::comm_world(), properties);
    if (dbo == nullptr || !dbo->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    // NOTE: 'output_region' owns 'dbo' pointer at this time
    Ioss::Region output_region(dbo, "skin");
    output_region.property_add(Ioss::Property(std::string("code_name"), codename));
    output_region.property_add(Ioss::Property(std::string("code_version"), version));

    output_region.begin_mode(Ioss::STATE_DEFINE_MODEL);

    auto *nbo = new Ioss::NodeBlock(output_region.get_database(), "nodeblock_1", ref_count, 3);

    // Count number of nodes owned by this processor (owner_out[i] == myProcessor);
    size_t owned = std::count_if(owner_out.begin(), owner_out.end(),
                                 [my_rank](int i) { return i == my_rank; });
    nbo->property_add(Ioss::Property("locally_owned_count", (INT)owned));

    output_region.add(nbo);

    // Output element blocks: will have a "skin" block for each input element block.
    //
    // NOTE: currently only handle a single skin face topology per
    // input element block. Eventually may need to handle the case
    // where there are two skin face topologies per input element
    // block (wedge -> tri and quad).

    // Count faces per element block and create output element block...
    for (auto &eb : ebs) {
      const std::string &name      = eb->name();
      const auto        &boundary  = boundary_faces[name];
      auto               face_topo = eb->topology()->face_type(0);
      std::string        topo      = "shell";
      if (face_topo == nullptr) {
        std::ostringstream errmsg;
        fmt::print(errmsg,
                   "ERROR: Block '{}' with topology '{}' does not have"
                   " a unique face topology.\nThis is not supported at this time.\n",
                   name, eb->topology()->name());
        IOSS_ERROR(errmsg);
      }
      if (face_topo->name() == "tri3") {
        topo = "trishell";
      }
      auto block =
          new Ioss::ElementBlock(output_region.get_database(), name, topo, boundary.size());
      output_region.add(block);
    }
    output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

    output_region.begin_mode(Ioss::STATE_MODEL);
    nbo->put_field_data("ids", ref_ids);
    nbo->put_field_data("owning_processor", owner_out);
    nbo->put_field_data("mesh_model_coordinates", coord_out);
    Ioss::Utils::clear(coord_out);
    Ioss::Utils::clear(owner_out);
    Ioss::Utils::clear(ref_ids);

    bool use_face_hash_ids = interFace.useFaceHashIds_;
    INT  fid               = 0;
    for (auto &eb : ebs) {
      const std::string &name       = eb->name();
      const auto        &boundary   = boundary_faces[name];
      auto              *block      = output_region.get_element_block(name);
      size_t             node_count = block->topology()->number_corner_nodes();

      Ioss::IntVector conn;
      Ioss::IntVector elids;
      conn.reserve(node_count * boundary.size());
      elids.reserve(boundary.size());

      for (auto &face : boundary) {
        if (use_face_hash_ids) {
          fid = face.hashId_;
          if (fid < 0) { // Due to (unsigned)size_t -> (signed)INT conversion
            fid = -fid;
          }
        }
        else {
          fid = face.element[0];
        }

        for (size_t i = 0; i < node_count; i++) {
          conn.push_back(face.connectivity_[i]);
        }
        elids.push_back(fid);
      }
      block->put_field_data("ids", elids);
      block->put_field_data("connectivity", conn);
    }
    output_region.end_mode(Ioss::STATE_MODEL);

    size_t ts_count = region.get_optional_property("state_count", 0);
    if (ts_count > 0 && interFace.output_transient()) {

      // Transfer all nodal variable names to the output database.
      Ioss::NodeBlock *nbi = region.get_node_blocks()[0];

      output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
      Ioss::NameList fields = nbi->field_describe(Ioss::Field::TRANSIENT);
      for (const auto &field_name : fields) {
        Ioss::Field field = nbi->get_field(field_name);
        field.reset_count(ref_count);
        nbo->field_add(field);
      }

      // All element block variables...
      for (auto &eb : ebs) {
        const std::string &name = eb->name();
        auto              *ebo  = output_region.get_element_block(name);
        if (ebo != nullptr) {
          auto           ebo_count = ebo->entity_count();
          Ioss::NameList efields   = eb->field_describe(Ioss::Field::TRANSIENT);
          for (const auto &field_name : efields) {
            Ioss::Field field = eb->get_field(field_name);
            field.reset_count(ebo_count);
            ebo->field_add(field);
          }
        }
      }

      output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

      auto selected_steps = get_selected_steps(region, interFace);

      output_region.begin_mode(Ioss::STATE_TRANSIENT);
      for (size_t istep = 1; istep <= ts_count; istep++) {
        if (selected_steps[istep] != 1) {
          continue;
        }
        double time  = region.get_state_time(istep);
        int    ostep = output_region.add_state(time);
        fmt::print(stderr, "\r\tTime step {:5d} at time {:10.5e}", ostep, time);

        output_region.begin_state(ostep);
        region.begin_state(istep);

        transfer_field_data(nbi, nbo, Ioss::Field::TRANSIENT, ref_nodes);

        // All element block variables...
        for (auto &eb : ebs) {
          const std::string &name = eb->name();
          auto              *ebo  = output_region.get_element_block(name);
          if (ebo != nullptr) {
            transfer_field_data(eb, ebo, Ioss::Field::TRANSIENT, boundary_faces[name]);
          }
        }

        output_region.end_state(ostep);
        region.end_state(istep);
      }
      output_region.end_mode(Ioss::STATE_TRANSIENT);
    }

    if (my_rank == 0) {
      output_region.output_summary(std::cerr, false);
    }
  }

  Ioss::PropertyManager set_properties(Skinner::Interface &interFace)
  {
    Ioss::PropertyManager properties{};
    if (interFace.ints_64_bit()) {
      properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
      properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    }

    if (interFace.debug) {
      properties.add(Ioss::Property("LOGGING", 1));
    }

    if (!interFace.decomp_method.empty()) {
      properties.add(Ioss::Property("DECOMPOSITION_METHOD", interFace.decomp_method));
    }

    if (interFace.compression_level > 0 || interFace.shuffle) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
      properties.add(Ioss::Property("COMPRESSION_LEVEL", interFace.compression_level));
      properties.add(Ioss::Property("COMPRESSION_SHUFFLE", interFace.shuffle));
    }

    if (interFace.compose_output == "default") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "NO"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "NO"));
    }
    else if (interFace.compose_output == "external") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "NO"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "NO"));
    }
    else if (interFace.compose_output != "none") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "YES"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "YES"));
    }

    if (interFace.netcdf4_) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }

    return properties;
  }

  void transfer_field_data_internal(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                                    const std::string &field_name, const Ioss::IntVector &ref_nodes)
  {
    if (oge->field_exists(field_name)) {
      std::vector<double> in;
      ige->get_field_data(field_name, in);

      // Determine component count...
      auto                field      = ige->get_field(field_name);
      int                 comp_count = field.get_component_count(Ioss::Field::InOut::INPUT);
      std::vector<double> out;
      out.reserve(oge->entity_count() * comp_count);

      for (size_t i = 0; i < ref_nodes.size(); i++) {
        if (ref_nodes[i] == 1) {
          for (int j = 0; j < comp_count; j++) {
            auto value = in[comp_count * i + j];
            out.push_back(value);
          }
        }
      }
      assert(out.size() == (size_t)oge->entity_count() * (size_t)comp_count);
      oge->put_field_data(field_name, out);
    }
  }

  void transfer_field_data(Ioss::GroupingEntity *ige, Ioss::GroupingEntity *oge,
                           Ioss::Field::RoleType role, const Ioss::IntVector &ref_nodes)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields = ige->field_describe(role);

    for (const auto &field_name : state_fields) {
      if (oge->field_exists(field_name)) {
        transfer_field_data_internal(ige, oge, field_name, ref_nodes);
      }
    }
  }

  void transfer_field_data_internal(Ioss::EntityBlock *ieb, Ioss::EntityBlock *oeb,
                                    const std::string             &field_name,
                                    const std::vector<Ioss::Face> &boundary_faces)
  {
    if (oeb->field_exists(field_name)) {
      std::vector<double> in;
      ieb->get_field_data(field_name, in);

      // Determine component count...
      auto                field      = ieb->get_field(field_name);
      int                 comp_count = field.get_component_count(Ioss::Field::InOut::INPUT);
      std::vector<double> out;
      out.reserve(oeb->entity_count() * comp_count);

      auto offset = ieb->get_offset();
      for (const auto &face : boundary_faces) {
        auto element_id = face.element[0] / 10 - 1;
        for (int j = 0; j < comp_count; j++) {
          auto value = in[comp_count * (element_id - offset) + j];
          out.push_back(value);
        }
      }
      assert(out.size() == (size_t)oeb->entity_count() * (size_t)comp_count);
      oeb->put_field_data(field_name, out);
    }
  }

  void transfer_field_data(Ioss::EntityBlock *ieb, Ioss::EntityBlock *oeb,
                           Ioss::Field::RoleType          role,
                           const std::vector<Ioss::Face> &boundary_faces)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields = ieb->field_describe(role);

    for (const auto &field_name : state_fields) {
      if (oeb->field_exists(field_name)) {
        transfer_field_data_internal(ieb, oeb, field_name, boundary_faces);
      }
    }
  }

} // namespace
