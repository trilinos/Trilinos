// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "io_info.h"

#include "Ioss_Hex8.h"
#include "Ioss_Sort.h"
#include "Ioss_StructuredBlock.h"
#include "Ioss_ZoneConnectivity.h"
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <stdint.h>
#include <tokenize.h>

#define FMT_DEPRECATED_OSTREAM
#include <cstddef>
#include <cstdlib>
#include <fmt/format.h>
#include <iostream>
#include <string>
#include <vector>

#include "Ioss_Assembly.h"
#include "Ioss_Blob.h"
#include "Ioss_BoundingBox.h"
#include "Ioss_CoordinateFrame.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DataSize.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_EdgeBlock.h"
#include "Ioss_EdgeSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementSet.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_FaceBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_Field.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_IOFactory.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_SurfaceSplit.h"
#include "Ioss_Utils.h"
#include "Ioss_VariableType.h"
#if defined(SEACAS_HAVE_EXODUS)
#include "exodusII.h"
#endif
#include "info_interface.h"
#if defined(SEACAS_HAVE_CGNS)
#include <cgnslib.h>
#endif

// ========================================================================

namespace {
  void info_timesteps(Ioss::Region &region);
  void info_nodeblock(Ioss::Region &region, const Info::Interface &interFace);
  void info_edgeblock(Ioss::Region &region, const Info::Interface &interFace);
  void info_faceblock(Ioss::Region &region, const Info::Interface &interFace);
  void info_elementblock(Ioss::Region &region, const Info::Interface &interFace);
  void info_structuredblock(Ioss::Region &region, const Info::Interface &interFace);

  void info_nodesets(Ioss::Region &region, const Info::Interface &interFace);
  void info_edgesets(Ioss::Region &region, const Info::Interface &interFace);
  void info_facesets(Ioss::Region &region, const Info::Interface &interFace);
  void info_elementsets(Ioss::Region &region, const Info::Interface &interFace);

  void info_sidesets(Ioss::Region &region, const Info::Interface &interFace);
  void info_coordinate_frames(Ioss::Region &region);
  void info_assemblies(Ioss::Region &region);
  void info_region(Ioss::Region &region);
  void info_blobs(Ioss::Region &region);

  void info_aliases(const Ioss::Region &region, const Ioss::GroupingEntity *ige, bool nl_pre,
                    bool nl_post);

  void info_variable_types()
  {
    auto var_list = Ioss::VariableType::external_types(Ioss::VariableType::Type::UNKNOWN);
    for (auto &var : var_list) {
      var->print();
    }
  }

  void file_info(const Info::Interface &interFace);
  void group_info(Info::Interface &interFace);

  void info_df(const Ioss::GroupingEntity *ge, const std::string &prefix)
  {
    int64_t             num_dist = ge->get_property("distribution_factor_count").get_int();
    std::vector<double> df;
    // Do even if num_dist == 0 so parallel does not assert.
    if (ge->field_exists("distribution_factors")) {
      ge->get_field_data("distribution_factors", df);
    }
    if (num_dist > 0) {
      auto mm = std::minmax_element(df.begin(), df.end());
      fmt::print("{}Distribution Factors: ", prefix);
      if (*mm.first == *mm.second) {
        fmt::print("all values = {:#}\n", *mm.first);
      }
      else {
        fmt::print("minimum value = {:#}, maximum value = {:#}\n", *mm.first, *mm.second);
      }
    }
  }

  template <typename T> void print_bbox(const T &entity)
  {
    Ioss::AxisAlignedBoundingBox bbox = entity.get_bounding_box();
    fmt::print("\tBounding Box: Minimum X,Y,Z = {:12.4e}\t{:12.4e}\t{:12.4e}\n"
               "\t              Maximum X,Y,Z = {:12.4e}\t{:12.4e}\t{:12.4e}\n"
               "\t                Range X,Y,Z = {:12.4e}\t{:12.4e}\t{:12.4e}\n",
               bbox.xmin, bbox.ymin, bbox.zmin, bbox.xmax, bbox.ymax, bbox.zmax,
               bbox.xmax - bbox.xmin, bbox.ymax - bbox.ymin, bbox.zmax - bbox.zmin);
  }

  std::string name(const Ioss::GroupingEntity *entity)
  {
    return entity->type_string() + " '" + entity->name() + "'";
  }

  int64_t id(Ioss::GroupingEntity *entity)
  {
    int64_t id = entity->get_optional_property("id", -1);
    return id;
  }

  Ioss::PropertyManager set_properties(const Info::Interface &interFace)
  {
    Ioss::PropertyManager properties{};
    if (!interFace.decomp_method().empty()) {
      properties.add(Ioss::Property("DECOMPOSITION_METHOD", interFace.decomp_method()));
    }
    return properties;
  }
} // namespace

void hex_volume(Ioss::ElementBlock *block, const std::vector<double> &coordinates);

namespace {
  void element_volume(Ioss::Region &region)
  {
    std::vector<double> coordinates;
    Ioss::NodeBlock    *nb = region.get_node_blocks()[0];
    nb->get_field_data("mesh_model_coordinates", coordinates);

    const Ioss::ElementBlockContainer &ebs = region.get_element_blocks();
    for (auto &eb : ebs) {
      if (eb->topology()->name() == Ioss::Hex8::name) {
        hex_volume(eb, coordinates);
      }
    }
  }

#if defined(SEACAS_HAVE_EXODUS)
  int print_groups(int exoid, std::string prefix)
  {
    int   idum;
    float rdum;
    char  group_name[33];
    // Print name of this group...
    ex_inquire(exoid, EX_INQ_GROUP_NAME, &idum, &rdum, group_name);
    if (group_name[0] == '/') {
      fmt::print("{}/ (root)\n", prefix);
    }
    else {
      fmt::print("{}{}\n", prefix, group_name);
    }

    int              num_children = ex_inquire_int(exoid, EX_INQ_NUM_CHILD_GROUPS);
    std::vector<int> children(num_children);
    ex_get_group_ids(exoid, nullptr, Data(children));
    prefix += '\t';
    for (int i = 0; i < num_children; i++) {
      print_groups(children[i], prefix);
    }
    return 0;
  }
#endif

  void group_info(Info::Interface &interFace)
  {
#if defined(SEACAS_HAVE_EXODUS)
    // Assume exodusII...
    std::string inpfile       = interFace.filename();
    float       vers          = 0.0;
    int         CPU_word_size = 0;
    int         IO_word_size  = 0;

    int exoid = ex_open(inpfile.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &vers);

    print_groups(exoid, "\t");
#endif
  }

  void file_info(const Info::Interface &interFace)
  {
    std::string inpfile    = interFace.filename();
    std::string input_type = interFace.type();

    Ioss::PropertyManager properties = set_properties(interFace);

    const auto custom_field = interFace.custom_field();
    if (!custom_field.empty()) {
      auto suffices = Ioss::tokenize(custom_field, ",");
      if (suffices.size() > 1) {
        Ioss::VariableType::create_named_suffix_type("UserDefined", suffices);
      }
    }

    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    auto mode = interFace.query_timesteps_only() ? Ioss::QUERY_TIMESTEPS_ONLY : Ioss::READ_RESTART;
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, mode,
                                                    Ioss::ParallelUtils::comm_world(), properties);

    Ioss::io_info_set_db_properties(interFace, dbi);

    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");
    if (interFace.query_timesteps_only()) {
      info_timesteps(region);
    }
    else {
      Ioss::io_info_file_info(interFace, region);
    }
  }

  void info_nodeblock(Ioss::Region &region, const Ioss::NodeBlock &nb,
                      const Info::Interface &interFace, const std::string &prefix = "")
  {
    int64_t num_nodes  = nb.entity_count();
    int64_t num_attrib = nb.get_property("attribute_count").get_int();
    fmt::print("\n{}{} {:14} nodes, {:3d} attributes.\n", prefix, name(&nb),
               fmt::group_digits(num_nodes), num_attrib);
    if (interFace.check_node_status()) {
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
    Ioss::Utils::info_fields(&nb, Ioss::Field::MAP, prefix + "\tMap Fields: ", "\n\t\t" + prefix);
    Ioss::Utils::info_fields(&nb, Ioss::Field::ATTRIBUTE,
                             prefix + "\tAttributes: ", "\n\t\t" + prefix);
    Ioss::Utils::info_fields(&nb, Ioss::Field::TRANSIENT,
                             prefix + "\tTransient:  ", "\n\t\t" + prefix,
                             interFace.field_details());

    if (interFace.compute_bbox() && region.mesh_type() != Ioss::MeshType::STRUCTURED) {
      print_bbox(nb);
    }
  }

  void info_nodeblock(Ioss::Region &region, const Info::Interface &interFace)
  {
    const Ioss::NodeBlockContainer &nbs = region.get_node_blocks();
    for (const auto &nb : nbs) {
      info_nodeblock(region, *nb, interFace, "");
    }
  }

  void info_structuredblock(Ioss::Region &region, const Info::Interface &interFace)
  {
    bool                                  parallel = region.get_database()->is_parallel();
    const Ioss::StructuredBlockContainer &sbs      = region.get_structured_blocks();
    for (auto &sb : sbs) {
      int64_t num_cell = sb->get_property("cell_count").get_int();
      int64_t num_node = sb->get_property("node_count").get_int();

      fmt::print("\n{} {}", name(sb), fmt::join(sb->get_ijk_global(), "x"));

      if (parallel) {
        fmt::print(" [{}, Offset = {}] ", fmt::join(sb->get_ijk_local(), "x"),
                   fmt::join(sb->get_ijk_offset(), ", "));
      }

      fmt::print("  {:14} cells, {:14} nodes ", fmt::group_digits(num_cell),
                 fmt::group_digits(num_node));

      info_aliases(region, sb, true, false);
      Ioss::Utils::info_fields(sb, Ioss::Field::TRANSIENT, "\n\tTransient:  ", "\n\t",
                               interFace.field_details());
      Ioss::Utils::info_fields(sb, Ioss::Field::REDUCTION, "\n\tTransient (Reduction):  ", "\t");
      info_nodeblock(region, sb->get_node_block(), interFace, "\t");
      fmt::print("\n");

      if (!sb->m_zoneConnectivity.empty()) {
        fmt::print("\tConnectivity with other blocks:\n");
        for (const auto &zgc : sb->m_zoneConnectivity) {
          std::cout << zgc << "\n";
        }
      }
      if (!sb->m_boundaryConditions.empty()) {
        fmt::print("\tBoundary Conditions:\n");
        // NOTE: The sort here is just to make io_info more useful for regression testing.
        //       With the sort, we get more reproducible output.  For now, only needed for BC...
        auto sb_bc = sb->m_boundaryConditions;
        Ioss::sort(sb_bc.begin(), sb_bc.end(),
                   [](const Ioss::BoundaryCondition &a, const Ioss::BoundaryCondition &b) {
                     return a.m_bcName < b.m_bcName;
                   });

        for (const auto &bc : sb_bc) {
          std::cout << bc << "\n";
        }
      }
      if (interFace.compute_bbox()) {
        print_bbox(*sb);
      }
    }
  }

  void info_region(Ioss::Region &region)
  {
    fmt::print("\nRegion '{}' (global)\n", region.name());
    Ioss::Utils::info_property(&region, Ioss::Property::ATTRIBUTE,
                               "\tAttributes (Reduction): ", "\t");
    Ioss::Utils::info_fields(&region, Ioss::Field::REDUCTION, "\tTransient  (Reduction):  ", "\t");
    info_variable_types();
  }

  void info_assemblies(Ioss::Region &region)
  {
    const auto &assem = region.get_assemblies();
    for (auto &as : assem) {
      fmt::print("\n{} id: {:6d}, contains: {} member(s) of type {:>10s}.\n\tMembers: ", name(as),
                 id(as), as->member_count(), as->contains_string());
      for (const auto &mem : as->get_members()) {
        fmt::print("'{}', ", mem->name());
      }

      info_aliases(region, as, true, false);
      fmt::print("\n");
      Ioss::Utils::info_property(as, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ", "\t");
      Ioss::Utils::info_fields(as, Ioss::Field::REDUCTION, "\tTransient  (Reduction):  ", "\t");
    }
  }

  void info_blobs(Ioss::Region &region)
  {
    const auto &blobs = region.get_blobs();
    for (auto &blob : blobs) {
      fmt::print("\n{} id: {:6d}, contains: {} item(s).\n", name(blob), id(blob),
                 blob->entity_count());

      info_aliases(region, blob, true, false);
      fmt::print("\n");
      Ioss::Utils::info_property(blob, Ioss::Property::ATTRIBUTE,
                                 "\tAttributes (Reduction): ", "\t");
      Ioss::Utils::info_fields(blob, Ioss::Field::TRANSIENT, "\tTransient:  ", "\t");
      Ioss::Utils::info_fields(blob, Ioss::Field::REDUCTION, "\tTransient  (Reduction):  ", "\t");
    }
  }

  void info_elementblock(Ioss::Region &region, const Info::Interface &interFace)
  {
    const Ioss::ElementBlockContainer &ebs = region.get_element_blocks();
    for (auto &eb : ebs) {
      int64_t num_elem = eb->entity_count();

      std::string type       = eb->topology()->name();
      int64_t     num_attrib = eb->get_property("attribute_count").get_int();
      fmt::print("\n{} id: {:6d}, topology: {:>10s}, {:14} elements, {:3d} attributes.", name(eb),
                 id(eb), type, fmt::group_digits(num_elem), num_attrib);

      info_aliases(region, eb, true, false);
      fmt::print("\n");
      Ioss::Utils::info_fields(eb, Ioss::Field::MAP, "\n\tMap Fields: ");
      Ioss::Utils::info_fields(eb, Ioss::Field::ATTRIBUTE, "\n\tAttributes: ");
      Ioss::Utils::info_property(eb, Ioss::Property::ATTRIBUTE, "\tAttributes (Reduction): ", "\t");

      if (interFace.adjacencies()) {
        Ioss::NameList blocks = eb->get_block_adjacencies();
        fmt::print("\n\tAdjacent to  {} element block(s):\t", blocks.size());
        for (const auto &block : blocks) {
          fmt::print("{}  ", block);
        }
      }
      Ioss::Utils::info_fields(eb, Ioss::Field::TRANSIENT, "\n\tTransient:  ", "\n\t",
                               interFace.field_details());
      Ioss::Utils::info_fields(eb, Ioss::Field::REDUCTION, "\n\tTransient  (Reduction):  ");

      if (interFace.compute_bbox()) {
        print_bbox(*eb);
      }
    }
  }

  void info_edgeblock(Ioss::Region &region, const Info::Interface &interFace)
  {
    const Ioss::EdgeBlockContainer &ebs = region.get_edge_blocks();
    for (auto &eb : ebs) {
      int64_t num_edge = eb->entity_count();

      std::string type       = eb->topology()->name();
      int64_t     num_attrib = eb->get_property("attribute_count").get_int();
      fmt::print("\n{} id: {:6d}, topology: {:>10s}, {:14} edges, {:3d} attributes.\n", name(eb),
                 id(eb), type, fmt::group_digits(num_edge), num_attrib);

      info_aliases(region, eb, false, true);
      Ioss::Utils::info_fields(eb, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

#if 0
        Ioss::NameList blocks = eb->get_block_adjacencies();
        fmt::print("\tAdjacent to  {} edge block(s):\t", blocks.size());
        for (auto &block : blocks) {
          fmt::print("{}  ", block);
        }
#endif
      Ioss::Utils::info_fields(eb, Ioss::Field::TRANSIENT, "\n\tTransient:  ", "\n\t",
                               interFace.field_details());
      Ioss::Utils::info_fields(eb, Ioss::Field::REDUCTION, "\n\tTransient (Reduction):  ");
      fmt::print("\n");
    }
  }

  void info_faceblock(Ioss::Region &region, const Info::Interface &interFace)
  {
    const Ioss::FaceBlockContainer &ebs = region.get_face_blocks();
    for (auto &eb : ebs) {
      int64_t num_face = eb->entity_count();

      std::string type       = eb->topology()->name();
      int64_t     num_attrib = eb->get_property("attribute_count").get_int();
      fmt::print("\n{} id: {:6d}, topology: {:>10s}, {:14} faces, {:3d} attributes.\n", name(eb),
                 id(eb), type, fmt::group_digits(num_face), num_attrib);

      info_aliases(region, eb, false, true);
      Ioss::Utils::info_fields(eb, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

#if 0
        Ioss::NameList blocks = eb->get_block_adjacencies();
        fmt::print("\tAdjacent to  {} face block(s):\t", blocks.size());
        for (auto &block : blocks) {
          fmt::print("{}  ", block);
        }
#endif
      Ioss::Utils::info_fields(eb, Ioss::Field::TRANSIENT, "\n\tTransient:  ", "\n\t",
                               interFace.field_details());
      Ioss::Utils::info_fields(eb, Ioss::Field::REDUCTION, "\n\tTransient (Reduction):  ");
      fmt::print("\n");
    }
  }

  void info_sidesets(Ioss::Region &region, const Info::Interface &interFace)
  {
    const Ioss::SideSetContainer &fss = region.get_sidesets();
    for (auto &fs : fss) {
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
      Ioss::Utils::info_fields(fs, Ioss::Field::TRANSIENT, "\n\tTransient: ", "\n\t",
                               interFace.field_details());
      Ioss::Utils::info_fields(fs, Ioss::Field::REDUCTION, "\n\tTransient (Reduction):  ");
      if (interFace.adjacencies()) {
        Ioss::NameList blocks;
        fs->block_membership(blocks);
        fmt::print("\n\t\tTouches {} element block(s):\t", blocks.size());
        for (const auto &block : blocks) {
          fmt::print("{}  ", block);
        }
      }
      fmt::print("\n");
      const Ioss::SideBlockContainer &fbs = fs->get_side_blocks();
      for (auto &fb : fbs) {
        int64_t count      = fb->entity_count();
        int64_t num_attrib = fb->get_property("attribute_count").get_int();
        int64_t num_dist   = fb->get_property("distribution_factor_count").get_int();
        fmt::print("\t{}, {:8} sides, {:3d} attributes, {:8} distribution factors.\n", name(fb),
                   fmt::group_digits(count), num_attrib, fmt::group_digits(num_dist));
        info_df(fb, "\t\t");
        Ioss::Utils::info_fields(fb, Ioss::Field::TRANSIENT, "\t\tTransient: ", "\n\t\t",
                                 interFace.field_details());
        Ioss::Utils::info_fields(fb, Ioss::Field::REDUCTION,
                                 "\t\tTransient (Reduction):  ", "\n\t\t");
      }
    }
  }

  void info_nodesets(Ioss::Region &region, const Info::Interface &interFace)
  {
    const Ioss::NodeSetContainer &nss = region.get_nodesets();
    for (auto &ns : nss) {
      int64_t count      = ns->entity_count();
      int64_t num_attrib = ns->get_property("attribute_count").get_int();
      int64_t num_dist   = ns->get_property("distribution_factor_count").get_int();
      fmt::print("\n{} id: {:6d}, {:8} nodes, {:3d} attributes, {:8} distribution factors.\n",
                 name(ns), id(ns), fmt::group_digits(count), num_attrib,
                 fmt::group_digits(num_dist));
      info_aliases(region, ns, false, true);
      info_df(ns, "\t");
      Ioss::Utils::info_fields(ns, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
      Ioss::Utils::info_fields(ns, Ioss::Field::TRANSIENT, "\tTransient:  ", "\n\t",
                               interFace.field_details());
      Ioss::Utils::info_fields(ns, Ioss::Field::REDUCTION, "\tTransient (Reduction):  ");
    }
  }

  void info_edgesets(Ioss::Region &region, const Info::Interface &interFace)
  {
    const Ioss::EdgeSetContainer &nss = region.get_edgesets();
    for (auto &ns : nss) {
      int64_t count      = ns->entity_count();
      int64_t num_attrib = ns->get_property("attribute_count").get_int();
      fmt::print("\n{} id: {:6d}, {:8} edges, {:3d} attributes.\n", name(ns), id(ns),
                 fmt::group_digits(count), num_attrib);
      info_aliases(region, ns, false, true);
      info_df(ns, "\t");
      Ioss::Utils::info_fields(ns, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
      Ioss::Utils::info_fields(ns, Ioss::Field::TRANSIENT, "\tTransient:  ", "\n\t",
                               interFace.field_details());
      Ioss::Utils::info_fields(ns, Ioss::Field::REDUCTION, "\tTransient (Reduction):  ");
    }
  }

  void info_facesets(Ioss::Region &region, const Info::Interface &interFace)
  {
    const Ioss::FaceSetContainer &fss = region.get_facesets();
    for (auto &fs : fss) {
      int64_t count      = fs->entity_count();
      int64_t num_attrib = fs->get_property("attribute_count").get_int();
      fmt::print("\n{} id: {:6d}, {:8} faces, {:3d} attributes.\n", name(fs), id(fs),
                 fmt::group_digits(count), num_attrib);
      info_aliases(region, fs, false, true);
      info_df(fs, "\t");
      Ioss::Utils::info_fields(fs, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
      Ioss::Utils::info_fields(fs, Ioss::Field::TRANSIENT, "\tTransient:  ", "\n\t",
                               interFace.field_details());
      Ioss::Utils::info_fields(fs, Ioss::Field::REDUCTION, "\tTransient (Reduction):  ");
    }
  }

  void info_elementsets(Ioss::Region &region, const Info::Interface &interFace)
  {
    const Ioss::ElementSetContainer &ess = region.get_elementsets();
    for (auto &es : ess) {
      int64_t count = es->entity_count();
      fmt::print("\n{} id: {:6d}, {:8} elements.\n", name(es), id(es), fmt::group_digits(count));
      info_aliases(region, es, false, true);
      info_df(es, "\t");
      Ioss::Utils::info_fields(es, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
      Ioss::Utils::info_fields(es, Ioss::Field::TRANSIENT, "\tTransient:  ", "\n\t",
                               interFace.field_details());
      Ioss::Utils::info_fields(es, Ioss::Field::REDUCTION, "\tTransient (Reduction):  ");
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
    Ioss::NameList aliases;
    if (region.get_aliases(ige->name(), ige->type(), aliases) > 0) {
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

  void info_timesteps(Ioss::Region &region)
  {
    int                 step_count = (int)region.get_property("state_count").get_int();
    std::vector<double> steps(step_count);

    for (int step = 0; step < step_count; step++) {
      double db_time = region.get_state_time(step + 1);
      steps[step]    = db_time;
    }
    auto mm = std::minmax_element(steps.begin(), steps.end());

    fmt::print("\nThere are {} time steps on the database.\n", step_count);
    fmt::print("\tMinimum Time = {:12.6e}, Maximum Time = {:12.6e}\n\n", *mm.first, *mm.second);

    fmt::print("\tStep Times: {:12.6e}\n", fmt::join(steps, ", "));
  }

} // namespace

namespace Ioss {
  void io_info_file_info(const Info::Interface &interFace) { file_info(interFace); }
  void io_info_group_info(Info::Interface &interFace) { group_info(interFace); }

  void io_info_set_db_properties(const Info::Interface &interFace, Ioss::DatabaseIO *dbi)
  {
    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    if (interFace.use_generic_names()) {
      dbi->set_use_generic_canonical_name(true);
    }

    if (interFace.surface_split_scheme() != Ioss::SPLIT_INVALID) {
      dbi->set_surface_split_type(Ioss::int_to_surface_split(interFace.surface_split_scheme()));
    }
    dbi->set_field_separator(interFace.field_suffix_separator());
    dbi->set_field_recognition(!interFace.disable_field_recognition());
    if (interFace.ints_64_bit()) {
      dbi->set_int_byte_size_api(Ioss::USE_INT64_API);
    }

    if (!interFace.groupname().empty()) {
      bool success = dbi->open_group(interFace.groupname());
      if (!success) {
        std::string inpfile = interFace.filename();
        fmt::print("ERROR: Unable to open group '{}' in file '{}'\n", interFace.groupname(),
                   inpfile);
        return;
      }
    }
  }

  void io_info_file_info(const Info::Interface &interFace, Ioss::Region &region)
  {
    Ioss::DatabaseIO *dbi = region.get_database();

    if (dbi == nullptr || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    bool parallel      = region.get_database()->is_parallel();
    int  parallel_size = region.get_database()->parallel_size();
    int  parallel_rank = region.get_database()->parallel_rank();
    region.get_database()->set_parallel_consistency(false);

    for (int proc = 0; proc < parallel_size; proc++) {
      if (proc == parallel_rank) {
        if (parallel) {
          fmt::print("\nProcessor {}", proc);
        }
        region.output_summary(std::cout, true);

        if (!interFace.summary()) {

          info_region(region);
          info_assemblies(region);
          info_nodeblock(region, interFace);
          info_edgeblock(region, interFace);
          info_faceblock(region, interFace);
          info_elementblock(region, interFace);
          info_structuredblock(region, interFace);

          info_nodesets(region, interFace);
          info_edgesets(region, interFace);
          info_facesets(region, interFace);
          info_elementsets(region, interFace);

          info_sidesets(region, interFace);
          info_blobs(region);
          info_coordinate_frames(region);
        }
      }
      region.get_database()->util().barrier();
    }

    if (interFace.compute_volume()) {
      element_volume(region);
    }
  }
} // namespace Ioss
