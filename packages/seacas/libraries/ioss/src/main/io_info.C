// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "io_info.h"
#include <Ioss_Hex8.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#if defined(SEACAS_HAVE_CGNS)
#include <cgnslib.h>
#endif

// ========================================================================

namespace {

  // Data space shared by most field input/output routines...
  std::vector<char> data;

  void info_nodeblock(Ioss::Region &region, const Info::Interface &interface);
  void info_edgeblock(Ioss::Region &region);
  void info_faceblock(Ioss::Region &region);
  void info_elementblock(Ioss::Region &region, const Info::Interface &interface);
  void info_structuredblock(Ioss::Region &region, const Info::Interface &interface);

  void info_nodesets(Ioss::Region &region);
  void info_edgesets(Ioss::Region &region);
  void info_facesets(Ioss::Region &region);
  void info_elementsets(Ioss::Region &region);

  void info_sidesets(Ioss::Region &region, const Info::Interface &interface);
  void info_coordinate_frames(Ioss::Region &region);

  void info_aliases(const Ioss::Region &region, const Ioss::GroupingEntity *ige, bool nl_pre,
                    bool nl_post);

  void info_fields(const Ioss::GroupingEntity *ige, Ioss::Field::RoleType role,
                   const std::string &header);

  void file_info(const Info::Interface &interface);
  void group_info(Info::Interface &interface);

  void info_df(const Ioss::GroupingEntity *ge, const std::string &prefix)
  {
    int64_t num_dist = ge->get_property("distribution_factor_count").get_int();
    if (num_dist > 0) {
      std::vector<double> df;
      ge->get_field_data("distribution_factors", df);
      auto mm = std::minmax_element(df.begin(), df.end());
      fmt::print("{}Distribution Factors: ", prefix);
      if (*mm.first == *mm.second) {
        fmt::print("all values = {}\n", *mm.first);
      }
      else {
        fmt::print("minimum value = {}, maximum value = {}\n", *mm.first, *mm.second);
      }
    }
  }

  std::string name(const Ioss::GroupingEntity *entity)
  {
    return entity->type_string() + " '" + entity->name() + "'";
  }

  int64_t id(Ioss::GroupingEntity *entity)
  {
    int64_t id = -1;
    if (entity->property_exists("id")) {
      id = entity->get_property("id").get_int();
    }
    return id;
  }

  Ioss::PropertyManager set_properties(const Info::Interface &interface)
  {
    Ioss::PropertyManager properties;
    if (!interface.decomp_method().empty()) {
      properties.add(Ioss::Property("DECOMPOSITION_METHOD", interface.decomp_method()));
    }
    return properties;
  }
} // namespace

void hex_volume(Ioss::ElementBlock *block, const std::vector<double> &coordinates);

namespace {
  void element_volume(Ioss::Region &region)
  {
    std::vector<double> coordinates;
    Ioss::NodeBlock *   nb = region.get_node_blocks()[0];
    nb->get_field_data("mesh_model_coordinates", coordinates);

    const Ioss::ElementBlockContainer &ebs = region.get_element_blocks();
    for (auto eb : ebs) {
      if (eb->get_property("topology_type").get_string() == Ioss::Hex8::name) {
        hex_volume(eb, coordinates);
      }
    }
  }

  int print_groups(int exoid, std::string prefix)
  {
#if defined(SEACAS_HAVE_EXODUS)
    int   idum;
    float rdum;
    char  group_name[33];
    // Print name of this group...
    ex_inquire(exoid, EX_INQ_GROUP_NAME, &idum, &rdum, group_name);
    fmt::print("{}{}\n", prefix, group_name);

    int              num_children = ex_inquire_int(exoid, EX_INQ_NUM_CHILD_GROUPS);
    std::vector<int> children(num_children);
    ex_get_group_ids(exoid, nullptr, TOPTR(children));
    prefix += '\t';
    for (int i = 0; i < num_children; i++) {
      print_groups(children[i], prefix);
    }
#endif
    return 0;
  }

  void group_info(Info::Interface &interface)
  {
#if defined(SEACAS_HAVE_EXODUS)
    // Assume exodusII...
    std::string inpfile       = interface.filename();
    float       vers          = 0.0;
    int         CPU_word_size = 0;
    int         IO_word_size  = 0;

    int exoid = ex_open(inpfile.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &vers);

    print_groups(exoid, "");
#endif
  }

  void file_info(const Info::Interface &interface)
  {
    std::string inpfile    = interface.filename();
    std::string input_type = interface.type();

    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::PropertyManager properties = set_properties(interface);

    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART,
                                                    (MPI_Comm)MPI_COMM_WORLD, properties);

    Ioss::io_info_set_db_properties(interface, dbi);

    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");
    Ioss::io_info_file_info(interface, region);
  }

  void info_nodeblock(Ioss::Region &region, const Ioss::NodeBlock &nb,
                      const Info::Interface &interface, const std::string &prefix = "")
  {
    int64_t num_nodes  = nb.entity_count();
    int64_t num_attrib = nb.get_property("attribute_count").get_int();
    fmt::print("\n{}{} {:14n} nodes, {:3d} attributes.\n", prefix, name(&nb), num_nodes,
               num_attrib);
    if (interface.check_node_status()) {
      std::vector<char>    node_status;
      std::vector<int64_t> ids;
      nb.get_field_data("node_connectivity_status", node_status);
      nb.get_field_data("ids", ids);
      bool header = false;
      for (size_t j = 0; j < node_status.size(); j++) {
        if (node_status[j] == 0) {
          if (!header) {
            header = true;
            fmt::print("\t{}Unconnected nodes: {}", prefix, ids[j]);
          }
          else {
            fmt::print(", {}", ids[j]);
          }
        }
      }
      if (header) {
        fmt::print("\n");
      }
    }
    if (!nb.is_nonglobal_nodeblock()) {
      info_aliases(region, &nb, false, true);
    }
    info_fields(&nb, Ioss::Field::ATTRIBUTE, prefix + "\tAttributes: ");
    info_fields(&nb, Ioss::Field::TRANSIENT, prefix + "\tTransient: ");
  }

  void info_nodeblock(Ioss::Region &region, const Info::Interface &interface)
  {
    const Ioss::NodeBlockContainer &nbs = region.get_node_blocks();
    for (auto nb : nbs) {
      info_nodeblock(region, *nb, interface, "");
    }
  }

  void info_structuredblock(Ioss::Region &region, const Info::Interface &interface)
  {
    bool parallel      = region.get_database()->is_parallel();
    int  parallel_size = region.get_database()->parallel_size();

    const Ioss::StructuredBlockContainer &sbs = region.get_structured_blocks();
    for (int proc = 0; proc < parallel_size; proc++) {
      if (proc == region.get_database()->parallel_rank()) {
        if (parallel)
          fmt::print("\nProcessor {}", proc);
      }
      for (auto sb : sbs) {
        int64_t num_cell = sb->get_property("cell_count").get_int();
        int64_t num_node = sb->get_property("node_count").get_int();
        int64_t num_dim  = sb->get_property("component_degree").get_int();

        fmt::print("\n{} {}", name(sb), sb->get_property("ni_global").get_int());
        if (num_dim > 1) {
          fmt::print("x{}", sb->get_property("nj_global").get_int());
        }
        if (num_dim > 2) {
          fmt::print("x{}", sb->get_property("nk_global").get_int());
        }

        if (parallel) {
          fmt::print(" [{}x{}x{}, Offset = {}, {}, {}] ", sb->get_property("ni").get_int(),
                     sb->get_property("nj").get_int(), sb->get_property("nk").get_int(),
                     sb->get_property("offset_i").get_int(), sb->get_property("offset_j").get_int(),
                     sb->get_property("offset_k").get_int());
        }

        fmt::print("{:14n} cells, {:14n} nodes ", num_cell, num_node);

        info_aliases(region, sb, true, false);
        info_fields(sb, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
        info_nodeblock(region, sb->get_node_block(), interface, "\t");
        fmt::print("\n");

        if (!sb->m_zoneConnectivity.empty()) {
          fmt::print("\tConnectivity with other blocks:\n");
          for (const auto &zgc : sb->m_zoneConnectivity) {
            fmt::print("{}\n", zgc);
          }
        }
        if (!sb->m_boundaryConditions.empty()) {
          fmt::print("\tBoundary Conditions:\n");
          for (const auto &bc : sb->m_boundaryConditions) {
            fmt::print("{}\n", bc);
          }
        }
        if (interface.compute_bbox()) {
          Ioss::AxisAlignedBoundingBox bbox = sb->get_bounding_box();
          fmt::print("\tBounding Box: Minimum X,Y,Z = {:12.4e}\t{:12.4e}\t{:12.4e}\n"
                     "\t              Maximum X,Y,Z = {:12.4e}\t{:12.4e}\t{:12.4e}\n",
                     bbox.xmin, bbox.ymin, bbox.zmin, bbox.xmax, bbox.ymax, bbox.zmax);
        }
      }
    }
  }

  void info_elementblock(Ioss::Region &region, const Info::Interface &interface)
  {
    const Ioss::ElementBlockContainer &ebs = region.get_element_blocks();
    for (auto eb : ebs) {
      int64_t num_elem = eb->entity_count();

      std::string type       = eb->get_property("topology_type").get_string();
      int64_t     num_attrib = eb->get_property("attribute_count").get_int();
      fmt::print("\n{} id: {:6d}, topology: {:>10s}, {:14n} elements, {:3d} attributes.", name(eb),
                 id(eb), type, num_elem, num_attrib);

      info_aliases(region, eb, true, false);
      info_fields(eb, Ioss::Field::ATTRIBUTE, "\n\tAttributes: ");

      if (interface.adjacencies()) {
        std::vector<std::string> blocks;
        eb->get_block_adjacencies(blocks);
        fmt::print("\n\tAdjacent to  {} element block(s):\t", blocks.size());
        for (const auto &block : blocks) {
          fmt::print("{}  ", block);
        }
      }
      info_fields(eb, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
      fmt::print("\n");

      if (interface.compute_bbox()) {
        Ioss::AxisAlignedBoundingBox bbox = eb->get_bounding_box();
        fmt::print("\tBounding Box: Minimum X,Y,Z = {:12.4e}\t{:12.4e}\t{:12.4e}\n"
                   "\t              Maximum X,Y,Z = {:12.4e}\t{:12.4e}\t{:12.4e}\n",
                   bbox.xmin, bbox.ymin, bbox.zmin, bbox.xmax, bbox.ymax, bbox.zmax);
      }
    }
  }

  void info_edgeblock(Ioss::Region &region)
  {
    const Ioss::EdgeBlockContainer &ebs = region.get_edge_blocks();
    for (auto eb : ebs) {
      int64_t num_edge = eb->entity_count();

      std::string type       = eb->get_property("topology_type").get_string();
      int64_t     num_attrib = eb->get_property("attribute_count").get_int();
      fmt::print("\n{} id: {:6d}, topology: {:>10s}, {:14n} edges, {:3d} attributes.\n", name(eb),
                 id(eb), type, num_edge, num_attrib);

      info_aliases(region, eb, false, true);
      info_fields(eb, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

#if 0
        std::vector<std::string> blocks;
        eb->get_block_adjacencies(blocks);
        fmt::print("\tAdjacent to  {} edge block(s):\t", blocks.size());
        for (auto block : blocks) {
          fmt::print("{}  ", block);
        }
#endif
      info_fields(eb, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
      fmt::print("\n");
    }
  }

  void info_faceblock(Ioss::Region &region)
  {
    const Ioss::FaceBlockContainer &ebs = region.get_face_blocks();
    for (auto eb : ebs) {
      int64_t num_face = eb->entity_count();

      std::string type       = eb->get_property("topology_type").get_string();
      int64_t     num_attrib = eb->get_property("attribute_count").get_int();
      fmt::print("\n{} id: {:6d}, topology: {:>10s}, {:14n} faces, {:3d} attributes.\n", name(eb),
                 id(eb), type, num_face, num_attrib);

      info_aliases(region, eb, false, true);
      info_fields(eb, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

#if 0
        std::vector<std::string> blocks;
        eb->get_block_adjacencies(blocks);
        fmt::print("\tAdjacent to  {} face block(s):\t", blocks.size());
        for (auto block : blocks) {
          fmt::print("{}  ", block);
        }
#endif
      info_fields(eb, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
      fmt::print("\n");
    }
  }

  void info_sidesets(Ioss::Region &region, const Info::Interface &interface)
  {
    const Ioss::SideSetContainer &fss = region.get_sidesets();
    for (auto fs : fss) {
      fmt::print("\n{} id: {:6d}", name(fs), id(fs));
      if (fs->property_exists("bc_type")) {
#if defined(SEACAS_HAVE_CGNS)
        auto bc_type = fs->get_property("bc_type").get_int();
        fmt::print(", boundary condition type: {} ({})", BCTypeName[bc_type], bc_type);
#else
        fmt::print(", boundary condition type: {}", fs->get_property("bc_type").get_int());
#endif
      }
      info_aliases(region, fs, true, false);
      info_fields(fs, Ioss::Field::TRANSIENT, "\n\tTransient: ");
      if (interface.adjacencies()) {
        std::vector<std::string> blocks;
        fs->block_membership(blocks);
        fmt::print("\n\tTouches {} element block(s):\t", blocks.size());
        for (const auto &block : blocks) {
          fmt::print("{}  ", block);
        }
        fmt::print("\n");
      }
      const Ioss::SideBlockContainer &fbs = fs->get_side_blocks();
      for (auto fb : fbs) {
        fmt::print("\n\t{}", name(fb));
        info_fields(fb, Ioss::Field::TRANSIENT, "\n\t\tTransient: ");
      }
    }
  }

  void info_nodesets(Ioss::Region &region)
  {
    const Ioss::NodeSetContainer &nss = region.get_nodesets();
    for (auto ns : nss) {
      int64_t count      = ns->entity_count();
      int64_t num_attrib = ns->get_property("attribute_count").get_int();
      int64_t num_dist   = ns->get_property("distribution_factor_count").get_int();
      fmt::print("\n{} id: {:6d}, {:8n} nodes, {:3d} attributes, {:8n} distribution factors.\n",
                 name(ns), id(ns), count, num_attrib, num_dist);
      info_aliases(region, ns, false, true);
      info_df(ns, "\t");
      info_fields(ns, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
      info_fields(ns, Ioss::Field::TRANSIENT, "\tTransient:  ");
    }
  }

  void info_edgesets(Ioss::Region &region)
  {
    const Ioss::EdgeSetContainer &nss = region.get_edgesets();
    for (auto ns : nss) {
      int64_t count      = ns->entity_count();
      int64_t num_attrib = ns->get_property("attribute_count").get_int();
      fmt::print("\n{} id: {:6d}, {:8n} edges, {:3d} attributes.\n", name(ns), id(ns), count,
                 num_attrib);
      info_aliases(region, ns, false, true);
      info_df(ns, "\t");
      info_fields(ns, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
      info_fields(ns, Ioss::Field::TRANSIENT, "\tTransient:  ");
    }
  }

  void info_facesets(Ioss::Region &region)
  {
    const Ioss::FaceSetContainer &fss = region.get_facesets();
    for (auto fs : fss) {
      int64_t count      = fs->entity_count();
      int64_t num_attrib = fs->get_property("attribute_count").get_int();
      fmt::print("\n{} id: {:6d}, {:8n} faces, {:3d} attributes.\n", name(fs), id(fs), count,
                 num_attrib);
      info_aliases(region, fs, false, true);
      info_df(fs, "\t");
      info_fields(fs, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
      info_fields(fs, Ioss::Field::TRANSIENT, "\tTransient:  ");
    }
  }

  void info_elementsets(Ioss::Region &region)
  {
    const Ioss::ElementSetContainer &ess = region.get_elementsets();
    for (auto es : ess) {
      int64_t count = es->entity_count();
      fmt::print("\n{} id: {:6d}, {:8n} elements.\n", name(es), id(es), count);
      info_aliases(region, es, false, true);
      info_df(es, "\t");
      info_fields(es, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
      info_fields(es, Ioss::Field::TRANSIENT, "\tTransient:  ");
    }
  }

  void info_coordinate_frames(Ioss::Region &region)
  {
    const Ioss::CoordinateFrameContainer &cf = region.get_coordinate_frames();
    for (const auto &frame : cf) {
      const double *origin = frame.origin();
      const double *a3pt   = frame.axis_3_point();
      const double *p13pt  = frame.plane_1_3_point();

      fmt::print("\nCoordinate Frame id: {:6d}, type tag '{}'\n"
                 "\tOrigin:          {}\t{}\t{}\n"
                 "\tAxis 3 Point:    {}\t{}\t{}\n"
                 "\tPlane 1-3 Point: {}\t{}\t{}\n",
                 frame.id(), frame.tag(), origin[0], origin[1], origin[2], a3pt[0], a3pt[1],
                 a3pt[2], p13pt[0], p13pt[1], p13pt[2]);
    }
  }

  void info_aliases(const Ioss::Region &region, const Ioss::GroupingEntity *ige, bool nl_pre,
                    bool nl_post)
  {
    std::vector<std::string> aliases;
    if (region.get_aliases(ige->name(), aliases) > 0) {
      if (nl_pre) {
        fmt::print("\n");
      }
      fmt::print("\tAliases: ");
      for (size_t i = 0; i < aliases.size(); i++) {
        if (aliases[i] != ige->name()) {
          if (i > 0) {
            fmt::print(", ");
          }
          fmt::print("{}", aliases[i]);
        }
      }
      if (nl_post) {
        fmt::print("\n");
      }
    }
  }

  void info_fields(const Ioss::GroupingEntity *ige, Ioss::Field::RoleType role,
                   const std::string &header)
  {
    Ioss::NameList fields;
    ige->field_describe(role, &fields);

    if (fields.empty()) {
      return;
    }

    if (!header.empty()) {
      fmt::print("{}", header);
    }
    // Iterate through results fields and transfer to output
    // database...
    for (const auto &field_name : fields) {
      const Ioss::VariableType *var_type   = ige->get_field(field_name).raw_storage();
      int                       comp_count = var_type->component_count();
      fmt::print("{:>16s}:{} ", field_name, comp_count);
    }
    if (!header.empty()) {
      fmt::print("\n");
    }
  }
} // namespace

namespace Ioss {
  void io_info_file_info(const Info::Interface &interface) { file_info(interface); }
  void io_info_group_info(Info::Interface &interface) { group_info(interface); }

  void io_info_set_db_properties(const Info::Interface &interface, Ioss::DatabaseIO *dbi)
  {
    std::string inpfile = interface.filename();

    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    if (interface.use_generic_names()) {
      dbi->set_use_generic_canonical_name(true);
    }

    dbi->set_surface_split_type(Ioss::int_to_surface_split(interface.surface_split_scheme()));
    dbi->set_field_separator(interface.field_suffix_separator());
    dbi->set_field_recognition(!interface.disable_field_recognition());
    if (interface.ints_64_bit()) {
      dbi->set_int_byte_size_api(Ioss::USE_INT64_API);
    }

    if (!interface.groupname().empty()) {
      bool success = dbi->open_group(interface.groupname());
      if (!success) {
        fmt::print("ERROR: Unable to open group '{}' in file '{}'\n", interface.groupname(),
                   inpfile);
        return;
      }
    }
  }

  void io_info_file_info(const Info::Interface &interface, Ioss::Region &region)
  {
    Ioss::DatabaseIO *dbi = region.get_database();

    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    // Get all properties of input database...
    region.output_summary(std::cerr, true);

    if (interface.summary() == 0) {
      info_nodeblock(region, interface);
      info_edgeblock(region);
      info_faceblock(region);
      info_elementblock(region, interface);
      info_structuredblock(region, interface);

      info_nodesets(region);
      info_edgesets(region);
      info_facesets(region);
      info_elementsets(region);

      info_sidesets(region, interface);
      info_coordinate_frames(region);
    }

    if (interface.compute_volume()) {
      element_volume(region);
    }
  }
} // namespace Ioss
