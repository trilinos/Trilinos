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
#if defined(SEACAS_HAVE_CGNS)
#include <cgnslib.h>
#endif

// ========================================================================

namespace {

  // Data space shared by most field input/output routines...
  std::vector<char> data;

  void info_nodeblock(Ioss::Region &region, const Info::Interface &interface, bool summary);
  void info_edgeblock(Ioss::Region &region, bool summary);
  void info_faceblock(Ioss::Region &region, bool summary);
  void info_elementblock(Ioss::Region &region, const Info::Interface &interface, bool summary);
  void info_structuredblock(Ioss::Region &region, const Info::Interface &interface, bool summary);

  void info_nodesets(Ioss::Region &region, bool summary);
  void info_edgesets(Ioss::Region &region, bool summary);
  void info_facesets(Ioss::Region &region, bool summary);
  void info_elementsets(Ioss::Region &region, bool summary);

  void info_sidesets(Ioss::Region &region, const Info::Interface &interface, bool summary);
  void info_commsets(Ioss::Region &region, bool summary);
  void info_coordinate_frames(Ioss::Region &region, bool summary);

  void info_aliases(Ioss::Region &region, Ioss::GroupingEntity *ige, bool nl_pre, bool nl_post);

  void info_fields(Ioss::GroupingEntity *ige, Ioss::Field::RoleType role,
                   const std::string &header);

  void info_properties(Ioss::GroupingEntity *ige);

  void file_info(const Info::Interface &interface);
  void group_info(Info::Interface &interface);

  std::string name(Ioss::GroupingEntity *entity)
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

    Ioss::ElementBlockContainer ebs = region.get_element_blocks();
    for (auto eb : ebs) {
      if (eb->get_property("topology_type").get_string() == Ioss::Hex8::name) {
        hex_volume(eb, coordinates);
      }
    }
  }

  int print_groups(int exoid, std::string prefix)
  {
#if !defined(NO_EXODUS_SUPPORT)
    int   idum;
    float rdum;
    char  group_name[33];
    // Print name of this group...
    ex_inquire(exoid, EX_INQ_GROUP_NAME, &idum, &rdum, group_name);
    OUTPUT << prefix << group_name << '\n';

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
#if !defined(NO_EXODUS_SUPPORT)
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

  void info_nodeblock(Ioss::Region &region, const Info::Interface &interface, bool summary)
  {
    Ioss::NodeBlockContainer nbs             = region.get_node_blocks();
    int64_t                  total_num_nodes = 0;
    if (summary) {
      int64_t degree = 0;
      for (auto nb : nbs) {
        int64_t num_nodes = nb->entity_count();
        total_num_nodes += num_nodes;
        degree = nb->get_property("component_degree").get_int();
      }
      OUTPUT << " Number of spatial dimensions =" << std::setw(12) << degree << "\n";
      OUTPUT << " Number of nodeblocks         =" << std::setw(12) << nbs.size() << "\t";
      OUTPUT << " Number of nodes            =" << std::setw(12) << total_num_nodes << "\n";
    }
    else {
      for (auto nb : nbs) {
        int64_t num_nodes  = nb->entity_count();
        int64_t num_attrib = nb->get_property("attribute_count").get_int();
        OUTPUT << '\n'
               << name(nb) << std::setw(12) << num_nodes << " nodes, " << std::setw(3) << num_attrib
               << " attributes.\n";
        if (interface.check_node_status()) {
          std::vector<char>    node_status;
          std::vector<int64_t> ids;
          nb->get_field_data("node_connectivity_status", node_status);
          nb->get_field_data("ids", ids);
          bool header = false;
          for (size_t j = 0; j < node_status.size(); j++) {
            if (node_status[j] == 0) {
              if (!header) {
                header = true;
                OUTPUT << "\tUnconnected nodes: " << ids[j];
              }
              else {
                OUTPUT << ", " << ids[j];
              }
            }
          }
          if (header) {
            OUTPUT << "\n";
          }
        }
        info_aliases(region, nb, false, true);
        info_fields(nb, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
        info_fields(nb, Ioss::Field::TRANSIENT, "\tTransient: ");
      }
    }
  }

  void info_structuredblock(Ioss::Region &region, const Info::Interface &interface, bool summary)
  {
    bool parallel      = region.get_database()->is_parallel();
    int  parallel_size = region.get_database()->parallel_size();

    int64_t total_cells = 0;
    int64_t total_nodes = 0;

    Ioss::StructuredBlockContainer sbs = region.get_structured_blocks();
    for (int proc = 0; proc < parallel_size; proc++) {
      if (proc == region.get_database()->parallel_rank()) {
        if (parallel && !summary) {
          OUTPUT << "\nProcessor " << proc;
        }
        for (auto sb : sbs) {
          int64_t num_cell = sb->get_property("cell_count").get_int();
          int64_t num_node = sb->get_property("node_count").get_int();
          int64_t num_dim  = sb->get_property("component_degree").get_int();

          total_cells += num_cell;
          total_nodes += num_node;

          if (!summary) {
            OUTPUT << '\n' << name(sb) << " " << sb->get_property("ni_global").get_int();
            if (num_dim > 1) {
              OUTPUT << "x" << sb->get_property("nj_global").get_int();
            }
            if (num_dim > 2) {
              OUTPUT << "x" << sb->get_property("nk_global").get_int();
            }

            if (parallel) {
              OUTPUT << " [" << sb->get_property("ni").get_int() << "x"
                     << sb->get_property("nj").get_int() << "x" << sb->get_property("nk").get_int()
                     << ", Offset = " << sb->get_property("offset_i").get_int() << ","
                     << sb->get_property("offset_j").get_int() << ","
                     << sb->get_property("offset_k").get_int() << "] ";
            }

            OUTPUT << std::setw(12) << num_cell << " cells, " << std::setw(12) << num_node
                   << " nodes ";

            info_aliases(region, sb, true, false);

            info_fields(sb, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
            OUTPUT << "\n";

            if (!sb->m_zoneConnectivity.empty()) {
              OUTPUT << "\tConnectivity with other blocks:\n";
              for (const auto &zgc : sb->m_zoneConnectivity) {
                OUTPUT << zgc << "\n";
              }
            }
            if (!sb->m_boundaryConditions.empty()) {
              OUTPUT << "\tBoundary Conditions:\n";
              for (const auto &bc : sb->m_boundaryConditions) {
                OUTPUT << bc << "\n";
              }
            }
            if (interface.compute_bbox()) {
              Ioss::AxisAlignedBoundingBox bbox = sb->get_bounding_box();
              OUTPUT << "\tBounding Box: Minimum X,Y,Z = " << std::setprecision(4)
                     << std::scientific << std::setw(12) << bbox.xmin << "\t" << std::setw(12)
                     << bbox.ymin << "\t" << std::setw(12) << bbox.zmin << "\n"
                     << "\t              Maximum X,Y,Z = " << std::setprecision(4)
                     << std::scientific << std::setw(12) << bbox.xmax << "\t" << std::setw(12)
                     << bbox.ymax << "\t" << std::setw(12) << bbox.zmax << "\n";
            }
          }
        }
      }
    }
    if (summary) {
      OUTPUT << " Number of structured blocks  =" << std::setw(12) << sbs.size() << "\t";
      OUTPUT << " Number of cells            =" << std::setw(12) << total_cells << "\n";
    }
  }

  void info_elementblock(Ioss::Region &region, const Info::Interface &interface, bool summary)
  {
    Ioss::ElementBlockContainer ebs            = region.get_element_blocks();
    int64_t                     total_elements = 0;
    for (auto eb : ebs) {
      int64_t num_elem = eb->entity_count();
      total_elements += num_elem;

      if (!summary) {
        std::string type       = eb->get_property("topology_type").get_string();
        int64_t     num_attrib = eb->get_property("attribute_count").get_int();
        OUTPUT << '\n'
               << name(eb) << " id: " << std::setw(6) << id(eb) << ", topology: " << std::setw(10)
               << type << ", " << std::setw(12) << num_elem << " elements, " << std::setw(3)
               << num_attrib << " attributes.";

        info_aliases(region, eb, true, false);
        info_fields(eb, Ioss::Field::ATTRIBUTE, "\n\tAttributes: ");

        if (interface.adjacencies()) {
          std::vector<std::string> blocks;
          eb->get_block_adjacencies(blocks);
          OUTPUT << "\n\tAdjacent to  " << blocks.size() << " element block(s):\t";
          for (const auto &block : blocks) {
            OUTPUT << block << "  ";
          }
        }
        info_fields(eb, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
        OUTPUT << "\n";

        if (interface.compute_bbox()) {
          Ioss::AxisAlignedBoundingBox bbox = eb->get_bounding_box();
          OUTPUT << "\tBounding Box: Minimum X,Y,Z = " << std::setprecision(4) << std::scientific
                 << std::setw(12) << bbox.xmin << "\t" << std::setw(12) << bbox.ymin << "\t"
                 << std::setw(12) << bbox.zmin << "\n"
                 << "\t              Maximum X,Y,Z = " << std::setprecision(4) << std::scientific
                 << std::setw(12) << bbox.xmax << "\t" << std::setw(12) << bbox.ymax << "\t"
                 << std::setw(12) << bbox.zmax << "\n";
        }
      }
    }
    if (summary) {
      OUTPUT << " Number of element blocks     =" << std::setw(12) << ebs.size() << "\t";
      OUTPUT << " Number of elements         =" << std::setw(12) << total_elements << "\n";
    }
  }

  void info_edgeblock(Ioss::Region &region, bool summary)
  {
    Ioss::EdgeBlockContainer ebs         = region.get_edge_blocks();
    int64_t                  total_edges = 0;
    for (auto eb : ebs) {
      int64_t num_edge = eb->entity_count();
      total_edges += num_edge;

      if (!summary) {
        std::string type       = eb->get_property("topology_type").get_string();
        int64_t     num_attrib = eb->get_property("attribute_count").get_int();
        OUTPUT << '\n'
               << name(eb) << " id: " << std::setw(6) << id(eb) << ", topology: " << std::setw(10)
               << type << ", " << std::setw(12) << num_edge << " edges, " << std::setw(3)
               << num_attrib << " attributes.\n";

        info_aliases(region, eb, false, true);
        info_fields(eb, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

#if 0
	std::vector<std::string> blocks;
	eb->get_block_adjacencies(blocks);
	OUTPUT << "\tAdjacent to  " << blocks.size() << " edge block(s):\t";
	for (auto block : blocks) {
	  OUTPUT << block << "  ";
	}
#endif
        info_fields(eb, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
        OUTPUT << "\n";
      }
    }
    if (summary) {
      OUTPUT << " Number of edge blocks        =" << std::setw(12) << ebs.size() << "\t";
      OUTPUT << " Number of edges            =" << std::setw(12) << total_edges << "\n";
    }
  }

  void info_faceblock(Ioss::Region &region, bool summary)
  {
    Ioss::FaceBlockContainer ebs         = region.get_face_blocks();
    int64_t                  total_faces = 0;
    for (auto eb : ebs) {
      int64_t num_face = eb->entity_count();
      total_faces += num_face;

      if (!summary) {
        std::string type       = eb->get_property("topology_type").get_string();
        int64_t     num_attrib = eb->get_property("attribute_count").get_int();
        OUTPUT << '\n'
               << name(eb) << " id: " << std::setw(6) << id(eb) << ", topology: " << std::setw(10)
               << type << ", " << std::setw(12) << num_face << " faces, " << std::setw(3)
               << num_attrib << " attributes.\n";

        info_aliases(region, eb, false, true);
        info_fields(eb, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

#if 0
	std::vector<std::string> blocks;
	eb->get_block_adjacencies(blocks);
	OUTPUT << "\tAdjacent to  " << blocks.size() << " face block(s):\t";
	for (auto block : blocks) {
	  OUTPUT << block << "  ";
	}
#endif
        info_fields(eb, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
        OUTPUT << "\n";
      }
    }
    if (summary) {
      OUTPUT << " Number of face blocks        =" << std::setw(12) << ebs.size() << "\t";
      OUTPUT << " Number of faces            =" << std::setw(12) << total_faces << "\n";
    }
  }

  void info_sidesets(Ioss::Region &region, const Info::Interface &interface, bool summary)
  {
    Ioss::SideSetContainer fss         = region.get_sidesets();
    int64_t                total_sides = 0;
    for (auto fs : fss) {
      if (!summary) {
        OUTPUT << '\n' << name(fs) << " id: " << std::setw(6) << id(fs);
        if (fs->property_exists("bc_type")) {
#if defined(SEACAS_HAVE_CGNS)
          auto bc_type = fs->get_property("bc_type").get_int();
          OUTPUT << ", boundary condition type: " << BCTypeName[bc_type] << " (" << bc_type << ")";
#else
          OUTPUT << ", boundary condition type: " << fs->get_property("bc_type").get_int();
#endif
        }
        info_aliases(region, fs, true, false);
        if (interface.adjacencies()) {
          std::vector<std::string> blocks;
          fs->block_membership(blocks);
          OUTPUT << "\n\tTouches " << blocks.size() << " element block(s):\t";
          for (const auto &block : blocks) {
            OUTPUT << block << "  ";
          }
          OUTPUT << "\n";
        }
      }
      if (!summary) {
        OUTPUT << "\n\tContains: \n";
      }

      Ioss::SideBlockContainer fbs = fs->get_side_blocks();
      for (auto fb : fbs) {
        int64_t num_side = fb->entity_count();
        if (!summary) {
          std::string fbtype  = fb->get_property("topology_type").get_string();
          std::string partype = fb->get_property("parent_topology_type").get_string();
          OUTPUT << "\t\t" << name(fb) << ", " << num_side << " " << fbtype << " sides"
                 << ", parent topology: " << partype;
          if (fb->parent_block() != nullptr) {
            const auto *parent = fb->parent_block();
            OUTPUT << ",\tparent block: '" << parent->name() << "' (" << parent->type_string()
                   << ")";
          }
          OUTPUT << "\n";
          if (interface.adjacencies()) {
            std::vector<std::string> blocks;
            fb->block_membership(blocks);
            OUTPUT << "\t\t\tTouches " << blocks.size() << " element block(s):\t";
            for (const auto &block : blocks) {
              OUTPUT << block << "  ";
            }
            OUTPUT << "\n";
          }
          OUTPUT << "\n";
          info_fields(fb, Ioss::Field::ATTRIBUTE, "\t\tAttributes: ");
          info_fields(fb, Ioss::Field::TRANSIENT, "\t\tTransient:  ");
        }
        total_sides += num_side;
      }
    }

    if (summary) {
      OUTPUT << " Number of side sets          =" << std::setw(12) << fss.size() << "\t";
      OUTPUT << " Number of element sides    =" << std::setw(12) << total_sides << "\n";
    }
  }

  void info_nodesets(Ioss::Region &region, bool summary)
  {
    Ioss::NodeSetContainer nss         = region.get_nodesets();
    int64_t                total_nodes = 0;
    for (auto ns : nss) {
      int64_t count      = ns->entity_count();
      int64_t num_attrib = ns->get_property("attribute_count").get_int();
      int64_t num_dist   = ns->get_property("distribution_factor_count").get_int();
      if (!summary) {
        OUTPUT << '\n'
               << name(ns) << " id: " << std::setw(6) << id(ns) << ", " << std::setw(8) << count
               << " nodes" << std::setw(3) << num_attrib << " attributes" << std::setw(8)
               << num_dist << " distribution factors.\n";
        info_aliases(region, ns, false, true);
        info_fields(ns, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
        info_fields(ns, Ioss::Field::TRANSIENT, "\tTransient:  ");
      }
      total_nodes += count;
    }
    if (summary) {
      OUTPUT << " Number of nodal point sets   =" << std::setw(12) << nss.size() << "\t";
      OUTPUT << " Length of node list        =" << std::setw(12) << total_nodes << "\n";
    }
  }

  void info_edgesets(Ioss::Region &region, bool summary)
  {
    Ioss::EdgeSetContainer nss         = region.get_edgesets();
    int64_t                total_edges = 0;
    for (auto ns : nss) {
      int64_t count      = ns->entity_count();
      int64_t num_attrib = ns->get_property("attribute_count").get_int();
      if (!summary) {
        OUTPUT << '\n'
               << name(ns) << " id: " << std::setw(6) << id(ns) << ", " << std::setw(8) << count
               << " edges" << std::setw(3) << num_attrib << " attributes.\n";
        info_aliases(region, ns, false, true);
        info_fields(ns, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
        info_fields(ns, Ioss::Field::TRANSIENT, "\tTransient:  ");
      }
      total_edges += count;
    }
    if (summary) {
      OUTPUT << " Number of edge sets          =" << std::setw(12) << nss.size() << "\t";
      OUTPUT << " Length of edge list        =" << std::setw(12) << total_edges << "\n";
    }
  }

  void info_facesets(Ioss::Region &region, bool summary)
  {
    Ioss::FaceSetContainer fss         = region.get_facesets();
    int64_t                total_faces = 0;
    for (auto fs : fss) {
      int64_t count      = fs->entity_count();
      int64_t num_attrib = fs->get_property("attribute_count").get_int();
      if (!summary) {
        OUTPUT << '\n'
               << name(fs) << " id: " << std::setw(6) << id(fs) << ", " << std::setw(8) << count
               << " faces" << std::setw(3) << num_attrib << " attributes.\n";
        info_aliases(region, fs, false, true);
        info_fields(fs, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
        info_fields(fs, Ioss::Field::TRANSIENT, "\tTransient:  ");
      }
      total_faces += count;
    }
    if (summary) {
      OUTPUT << " Number of face sets          =" << std::setw(12) << fss.size() << "\t";
      OUTPUT << " Length of face list        =" << std::setw(12) << total_faces << "\n";
    }
  }

  void info_elementsets(Ioss::Region &region, bool summary)
  {
    Ioss::ElementSetContainer ess            = region.get_elementsets();
    int64_t                   total_elements = 0;
    for (auto es : ess) {
      int64_t count = es->entity_count();
      if (!summary) {
        OUTPUT << '\n'
               << name(es) << " id: " << std::setw(6) << id(es) << ", " << std::setw(8) << count
               << " elements"
               << "\n";
        info_aliases(region, es, false, true);
        info_fields(es, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
        info_fields(es, Ioss::Field::TRANSIENT, "\tTransient:  ");
      }
      total_elements += count;
    }
    if (summary) {
      OUTPUT << " Number of element sets       =" << std::setw(12) << ess.size() << "\t";
      OUTPUT << " Length of element list     =" << std::setw(12) << total_elements << "\n";
    }
  }

  void info_commsets(Ioss::Region &region, bool summary)
  {
    // NOTE: This doesn't really do anything...
    Ioss::CommSetContainer css = region.get_commsets();
    for (auto cs : css) {
      std::string type = cs->get_property("entity_type").get_string();
    }
    OUTPUT << '\n';
  }

  void info_coordinate_frames(Ioss::Region &region, bool summary)
  {
    Ioss::CoordinateFrameContainer cf = region.get_coordinate_frames();
    for (const auto &frame : cf) {
      if (!summary) {
        const double *origin = frame.origin();
        const double *a3pt   = frame.axis_3_point();
        const double *p13pt  = frame.plane_1_3_point();

        OUTPUT << '\n'
               << "Coordinate Frame id: " << std::setw(6) << frame.id() << ", type tag '"
               << frame.tag() << "'\n"
               << "\tOrigin:          " << origin[0] << "\t" << origin[1] << "\t" << origin[2]
               << "\n"
               << "\tAxis 3 Point:    " << a3pt[0] << "\t" << a3pt[1] << "\t" << a3pt[2] << "\n"
               << "\tPlane 1-3 Point: " << p13pt[0] << "\t" << p13pt[1] << "\t" << p13pt[2] << "\n";
      }
    }
    if (summary) {
      OUTPUT << " Number of coordinate frames  =" << std::setw(12) << cf.size() << "\n";
    }
  }

  void info_aliases(Ioss::Region &region, Ioss::GroupingEntity *ige, bool nl_pre, bool nl_post)
  {
    std::vector<std::string> aliases;
    if (region.get_aliases(ige->name(), aliases) > 0) {
      if (nl_pre) {
        OUTPUT << "\n";
      }
      OUTPUT << "\tAliases: ";
      for (size_t i = 0; i < aliases.size(); i++) {
        if (aliases[i] != ige->name()) {
          if (i > 0) {
            OUTPUT << ", ";
          }
          OUTPUT << aliases[i];
        }
      }
      if (nl_post) {
        OUTPUT << "\n";
      }
    }
  }

  void info_fields(Ioss::GroupingEntity *ige, Ioss::Field::RoleType role, const std::string &header)
  {
    Ioss::NameList fields;
    ige->field_describe(role, &fields);

    if (fields.empty()) {
      return;
    }

    if (!header.empty()) {
      OUTPUT << header;
    }
    // Iterate through results fields and transfer to output
    // database...
    for (const auto &field_name : fields) {
      const Ioss::VariableType *var_type   = ige->get_field(field_name).raw_storage();
      int                       comp_count = var_type->component_count();
      OUTPUT << std::setw(16) << field_name << ":" << comp_count << " ";
    }
    if (!header.empty()) {
      OUTPUT << "\n";
    }
  }

  void info_properties(Ioss::GroupingEntity *ige)
  {
#if 0
    Ioss::NameList properties;
    ige->property_describe(&properties);

    // Iterate through properties and transfer to output database...
    for (auto property : properties) {
      OUTPUT << property << ", ";
    }
#endif
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
        OUTPUT << "ERROR: Unable to open group '" << interface.groupname() << "' in file '"
               << inpfile << "\n";
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
    bool summary = true;
    info_properties(&region);
    info_nodeblock(region, interface, summary);
    info_edgeblock(region, summary);
    info_faceblock(region, summary);
    info_elementblock(region, interface, summary);
    info_structuredblock(region, interface, summary);

    info_nodesets(region, summary);
    info_edgesets(region, summary);
    info_facesets(region, summary);
    info_elementsets(region, summary);

    info_sidesets(region, interface, summary);
    info_commsets(region, summary);
    info_coordinate_frames(region, summary);
    if (region.property_exists("state_count") && region.get_property("state_count").get_int() > 0) {
      std::pair<int, double> state_time_max = region.get_max_time();
      std::pair<int, double> state_time_min = region.get_min_time();
      OUTPUT << " Number of time steps on database     =" << std::setw(12)
             << region.get_property("state_count").get_int() << "\n"
             << "    Minimum time = " << state_time_min.second << " at step "
             << state_time_min.first << "\n"
             << "    Maximum time = " << state_time_max.second << " at step "
             << state_time_max.first << "\n\n";
    }

    if (interface.summary() == 0) {
      summary = false;
      info_properties(&region);
      info_nodeblock(region, interface, summary);
      info_edgeblock(region, summary);
      info_faceblock(region, summary);
      info_elementblock(region, interface, summary);
      info_structuredblock(region, interface, summary);

      info_nodesets(region, summary);
      info_edgesets(region, summary);
      info_facesets(region, summary);
      info_elementsets(region, summary);

      info_sidesets(region, interface, summary);
      info_commsets(region, summary);
      info_coordinate_frames(region, summary);
    }

    if (interface.compute_volume()) {
      element_volume(region);
    }
  }
} // namespace Ioss
