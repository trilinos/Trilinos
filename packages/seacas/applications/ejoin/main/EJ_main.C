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

namespace {

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

    unsigned int debug_level = interFace.debug();

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

namespace {
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
} // namespace
