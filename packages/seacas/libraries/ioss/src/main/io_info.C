// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
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

#include <Ioss_CodeTypes.h>

#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <string>
#include <cstring>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>
#include <Ioss_SurfaceSplit.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifndef NO_XDMF_SUPPORT
#include <xdmf/Ioxf_Initializer.h>
#endif

#define OUTPUT std::cerr

// ========================================================================

namespace {

  // Data space shared by most field input/output routines...
  std::vector<char> data;

  struct Globals
  {
    bool summary;
    bool check_node_status;
    double maximum_time;
    double minimum_time;
    int  surface_split_type;
    std::string working_directory;
  };

  void show_usage(const std::string &prog);
  void show_step(int istep, double time);

  void info_nodeblock(Ioss::Region &region, const Globals &options);
  void info_edgeblock(Ioss::Region &region, bool summary);
  void info_faceblock(Ioss::Region &region, bool summary);
  void info_elementblock(Ioss::Region &region, bool summary);

  void info_nodesets(Ioss::Region &region, bool summary);
  void info_edgesets(Ioss::Region &region, bool summary);
  void info_facesets(Ioss::Region &region, bool summary);
  void info_elementsets(Ioss::Region &region, bool summary);

  void info_sidesets(Ioss::Region &region, bool summary);
  void info_commsets(Ioss::Region &region, bool summary);

  void info_fields(Ioss::GroupingEntity *ige,
		   Ioss::Field::RoleType role,
		   const std::string &header);

  void info_field_data(Ioss::GroupingEntity *ige,
		       Ioss::Field::RoleType role,
		       bool do_connectivity = true);

  void info_properties(Ioss::GroupingEntity *ige);

  void info_field_data_internal(Ioss::GroupingEntity *ige,
				const std::string &field_name);

  void file_info(const std::string& inpfile, const std::string& input_type,
		 Globals& globals);

  std::string name(Ioss::GroupingEntity *entity) {
    return entity->type_string() + " " + entity->name();
  }

  int id(Ioss::GroupingEntity *entity) {
    int id = -1;
    if (entity->property_exists("id")) {
      id = entity->get_property("id").get_int();
    }
    return id;
  }


}
void hex_volume(Ioss::ElementBlock *block, std::vector<double> &coordinates);

// ========================================================================

namespace {
  std::string codename;
  std::string version = "1.0";
}

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  
  std::string in_type = "exodusII";

  Globals globals;

  globals.maximum_time = 0.0;
  globals.minimum_time = 0.0;
  globals.surface_split_type = 1;
  globals.check_node_status = false;
  
  codename = argv[0];
  size_t ind = codename.find_last_of("/", codename.size());
  if (ind != std::string::npos)
    codename = codename.substr(ind+1, codename.size());

  Ioss::Init::Initializer io;
#ifndef NO_XDMF_SUPPORT
  Ioxf::Initializer ioxf;
#endif

  // Skip past any options...
  int i=1;
  while (i < argc && argv[i][0] == '-') {
    if (std::strcmp("-directory", argv[i]) == 0 ||
	std::strcmp("-d", argv[i]) == 0) {
      i++;
      globals.working_directory = argv[i++];
    }
    else if (std::strcmp("--in_type", argv[i]) == 0) {
      i++;
      in_type = argv[i++];
    }
    else if (std::strcmp("-Maximum_Time", argv[i]) == 0) {
      i++;
      globals.maximum_time = std::strtod(argv[i++], NULL);
    }
    else if (std::strcmp("-Minimum_Time", argv[i]) == 0) {
      i++;
      globals.minimum_time = std::strtod(argv[i++], NULL);
    }
    else if (std::strcmp("-Surface_Split_Scheme", argv[i]) == 0) {
      i++;
      char *split_scheme = argv[i++];
      if (std::strcmp(split_scheme, "TOPOLOGY") == 0)
	globals.surface_split_type = 1;
      else if (std::strcmp(split_scheme, "ELEMENT_BLOCK") == 0)
	globals.surface_split_type = 2;
      else if (std::strcmp(split_scheme, "NO_SPLIT") == 0)
	globals.surface_split_type = 3;
    }

    else if (std::strcmp("-Node_Status", argv[i]) == 0) {
      i++;
      globals.check_node_status = true;
    }

    // Found an option.  See if it has an argument...
    else if (i+1 < argc && argv[i+1][0] == '-') {
      // No argument, another option
      i++;
    } else {
      // Skip the argument...
      i += 2;
    }
  }

  // Last argument is the filename...
  std::string in_file   = Ioss::Utils::local_filename(argv[argc-1], in_type, globals.working_directory);

  OUTPUT << "Input:    '" << in_file  << "', Type: " << in_type  << '\n';
  OUTPUT << '\n';

  file_info(in_file, in_type, globals);

  OUTPUT << "\n" << codename << " execution successful.\n";
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}

namespace {
  void show_usage(const std::string &prog)
  {
    OUTPUT << "\nUSAGE: " << prog << " input_database\n";
    OUTPUT << "       version: " << version << "\n";
    Ioss::NameList db_types;
    Ioss::IOFactory::describe(&db_types);
    OUTPUT << "\nSupports database types:\n\t";
    for (Ioss::NameList::const_iterator IF = db_types.begin(); IF != db_types.end(); ++IF) {
      OUTPUT << *IF << "  ";
    }
    OUTPUT << "\n\n";
  }

  void element_volume(Ioss::Region &region)
  {
    std::vector<double> coordinates;
    Ioss::NodeBlock *nb = region.get_node_blocks()[0];
    nb->get_field_data("mesh_model_coordinates", coordinates);

    Ioss::ElementBlockContainer ebs = region.get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
    while (i != ebs.end()) {
      if ((*i)->get_property("topology_type").get_string() == "hex8") {
	hex_volume(*i, coordinates);
      }
      ++i;
    }
  }

  void file_info(const std::string& inpfile, const std::string& input_type, Globals& globals)
  {
    //========================================================================
    // INPUT ...
    // NOTE: The "READ_RESTART" mode ensures that the node and element ids will be mapped.
    //========================================================================
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(input_type, inpfile, Ioss::READ_RESTART,
						    (MPI_Comm)MPI_COMM_WORLD);
    if (dbi == NULL || !dbi->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    dbi->set_surface_split_type(Ioss::int_to_surface_split(globals.surface_split_type));
    //    dbi->set_field_separator(0);
    dbi->set_node_global_id_backward_compatibility(false);
    
    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");

    // Get all properties of input database...
    globals.summary = true;
    info_properties(&region);
    info_nodeblock(region,    globals);
    info_edgeblock(region,    globals.summary);
    info_faceblock(region,    globals.summary);
    info_elementblock(region, globals.summary);

    info_nodesets(region,     globals.summary);
    info_edgesets(region,     globals.summary);
    info_facesets(region,     globals.summary);
    info_elementsets(region,  globals.summary);

    info_sidesets(region,     globals.summary);
    info_commsets(region,     globals.summary);

    if (region.property_exists("state_count") && region.get_property("state_count").get_int() > 0) {
      std::pair<int, double> state_time_max = region.get_max_time();
      std::pair<int, double> state_time_min = region.get_min_time();
      OUTPUT << " Number of time steps on database     =" << std::setw(9)
	     << region.get_property("state_count").get_int() << "\n"
	     << "    Minimum time = " << state_time_min.second << " at step " << state_time_min.first << "\n"
	     << "    Maximum time = " << state_time_max.second << " at step " << state_time_max.first << "\n\n";
    }

    globals.summary = false;
    info_properties(&region);
    info_nodeblock(region,    globals);
    info_edgeblock(region,    globals.summary);
    info_faceblock(region,    globals.summary);
    info_elementblock(region, globals.summary);

    info_nodesets(region,     globals.summary);
    info_edgesets(region,     globals.summary);
    info_facesets(region,     globals.summary);
    info_elementsets(region,  globals.summary);

    info_sidesets(region,     globals.summary);
    info_commsets(region,     globals.summary);

    bool do_volume = true;
    if (do_volume) {
      element_volume(region);
    }
  }


  void info_nodeblock(Ioss::Region &region, const Globals &options)
  {
    Ioss::NodeBlockContainer    nbs = region.get_node_blocks();
    Ioss::NodeBlockContainer::const_iterator i = nbs.begin();
    while (i != nbs.end()) {
      //      std::string name      = (*i)->name();
      int    num_nodes = (*i)->get_property("entity_count").get_int();
      int    degree    = (*i)->get_property("component_degree").get_int();
      int    num_attrib= (*i)->get_property("attribute_count").get_int();
      if (options.summary) {
	OUTPUT << " Number of spatial dimensions         =" << std::setw(9) << degree << "\n";
	OUTPUT << " Number of nodeblocks                 =" << std::setw(9) << 1 << "\t";
	OUTPUT << " Number of nodes            =" << std::setw(9) << num_nodes << "\n";
      } else {
	OUTPUT << '\n' << name(*i) 
	       << std::setw(9) << num_nodes << " nodes, "
	       << std::setw(3) << num_attrib << " attributes.\n";
	if (options.check_node_status) {
	  std::vector<char> node_status;
	  std::vector<int>  ids;
	  (*i)->get_field_data("node_connectivity_status", node_status);
	  (*i)->get_field_data("ids", ids);
	  bool header = false;
	  for (size_t j=0; j < node_status.size(); j++) {
	    if (node_status[j] == 0) {
	      if (!header) {
		header = true;
		OUTPUT << "\tUnconnected nodes: " << ids[j];
	      } else {
		OUTPUT << ", " << ids[j];
	      }
	    }
	  }
	  if (header)
	    OUTPUT << "\n";
	}
	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
	info_fields(*i, Ioss::Field::TRANSIENT, "\tTransient: ");
      }
      ++i;
    }
  }

  void info_elementblock(Ioss::Region &region, bool summary)
  {
    Ioss::ElementBlockContainer ebs = region.get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
    int total_elements = 0;
    while (i != ebs.end()) {
      int    num_elem  = (*i)->get_property("entity_count").get_int();
      total_elements += num_elem;

      if (!summary) {
	std::string type      = (*i)->get_property("topology_type").get_string();
	int    num_attrib= (*i)->get_property("attribute_count").get_int();
	OUTPUT << '\n' << name(*i)
	       << " id: " << std::setw(6) << id(*i)
	       << ", topology: " << std::setw(10) << type << ", "
	       << std::setw(9) << num_elem << " elements, "
	       << std::setw(3) << num_attrib << " attributes.\n";

	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

	std::vector<std::string> blocks;
	(*i)->get_block_adjacencies(blocks);
	OUTPUT << "\tAdjacent to  " << blocks.size() << " element block(s):\t";
	std::vector<std::string>::iterator b = blocks.begin();
	while (b != blocks.end()) {
	  OUTPUT << *b++ << "  ";
	}
	info_fields(*i, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
	OUTPUT << "\n";
      }

      ++i;
    }
    if (summary) {
      OUTPUT << " Number of element blocks             =" << std::setw(9) << ebs.size() << "\t";
      OUTPUT << " Number of elements         =" << std::setw(9) << total_elements << "\n";
    }

  }

  void info_edgeblock(Ioss::Region &region, bool summary)
  {
    Ioss::EdgeBlockContainer ebs = region.get_edge_blocks();
    Ioss::EdgeBlockContainer::const_iterator i = ebs.begin();
    int total_edges = 0;
    while (i != ebs.end()) {
      int    num_edge  = (*i)->get_property("entity_count").get_int();
      total_edges += num_edge;

      if (!summary) {
	std::string type      = (*i)->get_property("topology_type").get_string();
	int    num_attrib= (*i)->get_property("attribute_count").get_int();
	OUTPUT << '\n' << name(*i)
	       << " id: " << std::setw(6) << id(*i)
	       << ", topology: " << std::setw(10) << type << ", "
	       << std::setw(9) << num_edge << " edges, "
	       << std::setw(3) << num_attrib << " attributes.\n";

	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

	std::vector<std::string> blocks;
	(*i)->get_block_adjacencies(blocks);
	OUTPUT << "\tAdjacent to  " << blocks.size() << " edge block(s):\t";
	std::vector<std::string>::iterator b = blocks.begin();
	while (b != blocks.end()) {
	  OUTPUT << *b++ << "  ";
	}
	info_fields(*i, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
	OUTPUT << "\n";
      }

      ++i;
    }
    if (summary) {
      OUTPUT << " Number of edge blocks                =" << std::setw(9) << ebs.size() << "\t";
      OUTPUT << " Number of edges            =" << std::setw(9) << total_edges << "\n";
    }

  }

  void info_faceblock(Ioss::Region &region, bool summary)
  {
    Ioss::FaceBlockContainer ebs = region.get_face_blocks();
    Ioss::FaceBlockContainer::const_iterator i = ebs.begin();
    int total_faces = 0;
    while (i != ebs.end()) {
      int    num_face  = (*i)->get_property("entity_count").get_int();
      total_faces += num_face;

      if (!summary) {
	std::string type      = (*i)->get_property("topology_type").get_string();
	int    num_attrib= (*i)->get_property("attribute_count").get_int();
	OUTPUT << '\n' << name(*i)
	       << " id: " << std::setw(6) << id(*i)
	       << ", topology: " << std::setw(10) << type << ", "
	       << std::setw(9) << num_face << " faces, "
	       << std::setw(3) << num_attrib << " attributes.\n";

	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

	std::vector<std::string> blocks;
	(*i)->get_block_adjacencies(blocks);
	OUTPUT << "\tAdjacent to  " << blocks.size() << " face block(s):\t";
	std::vector<std::string>::iterator b = blocks.begin();
	while (b != blocks.end()) {
	  OUTPUT << *b++ << "  ";
	}
	info_fields(*i, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
	OUTPUT << "\n";
      }

      ++i;
    }
    if (summary) {
      OUTPUT << " Number of face blocks                =" << std::setw(9) << ebs.size() << "\t";
      OUTPUT << " Number of faces            =" << std::setw(9) << total_faces << "\n";
    }

  }

  void info_sidesets(Ioss::Region &region, bool summary)
  {
    Ioss::SideSetContainer      fss = region.get_sidesets();
    Ioss::SideSetContainer::const_iterator i = fss.begin();
    int total_sides = 0;
    while (i != fss.end()) {
      if (!summary) {
	std::vector<std::string> blocks;
	(*i)->block_membership(blocks);
	OUTPUT << '\n' << (*i)->type_string() << " " << std::setw(16)  << (*i)->name()
	       << " id: " << std::setw(6)<< id(*i) << ":\n";
	OUTPUT << "\tTouches " << blocks.size() << " element block(s):\t";
	std::vector<std::string>::iterator b = blocks.begin();
	while (b != blocks.end()) {
	  OUTPUT << *b++ << "  ";
	}
	OUTPUT << "\n\tContains: \n";
      }

      Ioss::SideBlockContainer fbs = (*i)->get_side_blocks();
      Ioss::SideBlockContainer::const_iterator j = fbs.begin();
      while (j != fbs.end()) {
	int    num_side  = (*j)->get_property("entity_count").get_int();
	if (!summary) {
	  std::string fbtype    = (*j)->get_property("topology_type").get_string();
	  std::string partype   = (*j)->get_property("parent_topology_type").get_string();
	  OUTPUT << "\t\t"
		 << (*j)->type_string() << " " << (*j)->name()
		 << ", "
		 << num_side << " " << fbtype << " sides"
		 << ", parent topology: " << partype 
		 << "\n";
	  std::vector<std::string> blocks;
	  (*j)->block_membership(blocks);
	  OUTPUT << "\t\t\tTouches " << blocks.size() << " element block(s):\t";
	  std::vector<std::string>::iterator b = blocks.begin();
	  while (b != blocks.end()) {
	    OUTPUT << *b++ << "  ";
	  }
	  OUTPUT << "\n";
	  info_fields(*i, Ioss::Field::TRANSIENT, "\tTransient:  ");
	}
	total_sides += num_side;
	++j;
      }
      ++i;
    }
    if (summary) {
      OUTPUT << " Number of element side sets          =" << std::setw(9) << fss.size() << "\t";
      OUTPUT << " Number of element sides    =" << std::setw(9) << total_sides << "\n";
    }
  }

  void info_nodesets(Ioss::Region &region, bool summary)
  {
    Ioss::NodeSetContainer      nss = region.get_nodesets();
    Ioss::NodeSetContainer::const_iterator i = nss.begin();
    int total_nodes = 0;
    while (i != nss.end()) {
      int    count     = (*i)->get_property("entity_count").get_int();
      int    num_attrib= (*i)->get_property("attribute_count").get_int();
      if (!summary) {
	OUTPUT << '\n' << (*i)->type_string() << " " << std::setw(16)  << (*i)->name()
	       << " id: " << std::setw(6) << id(*i)   << ", "
	       << std::setw(8) << count << " nodes" 
	       << std::setw(3) << num_attrib << " attributes.\n";
	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
	info_fields(*i, Ioss::Field::TRANSIENT, "\tTransient:  ");
      }
      total_nodes += count;
      ++i;
    }
    if (summary) {
      OUTPUT << " Number of nodal point sets           =" << std::setw(9) << nss.size() << "\t";
      OUTPUT << " Length of node list        =" << std::setw(9) << total_nodes << "\n";
    }
  }

  void info_edgesets(Ioss::Region &region, bool summary)
  {
    Ioss::EdgeSetContainer      nss = region.get_edgesets();
    Ioss::EdgeSetContainer::const_iterator i = nss.begin();
    int total_edges = 0;
    while (i != nss.end()) {
      int    count     = (*i)->get_property("entity_count").get_int();
      int    num_attrib= (*i)->get_property("attribute_count").get_int();
      if (!summary) {
	OUTPUT << '\n' << (*i)->type_string() << " " << std::setw(16)  << (*i)->name()
	       << " id: " << std::setw(6) << id(*i)   << ", "
	       << std::setw(8) << count << " edges"
	       << std::setw(3) << num_attrib << " attributes.\n";
	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
	info_fields(*i, Ioss::Field::TRANSIENT, "\tTransient:  ");
      }
      total_edges += count;
      ++i;
    }
    if (summary) {
      OUTPUT << " Number of edge sets                  =" << std::setw(9) << nss.size() << "\t";
      OUTPUT << " Length of edge list        =" << std::setw(9) << total_edges << "\n";
    }
  }

  void info_facesets(Ioss::Region &region, bool summary)
  {
    Ioss::FaceSetContainer      nss = region.get_facesets();
    Ioss::FaceSetContainer::const_iterator i = nss.begin();
    int total_faces = 0;
    while (i != nss.end()) {
      int    count     = (*i)->get_property("entity_count").get_int();
      int    num_attrib= (*i)->get_property("attribute_count").get_int();
      if (!summary) {
	OUTPUT << '\n' << (*i)->type_string() << " " << std::setw(16)  << (*i)->name()
	       << " id: " << std::setw(6) << id(*i)   << ", "
	       << std::setw(8) << count << " faces"
	       << std::setw(3) << num_attrib << " attributes.\n";
	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
	info_fields(*i, Ioss::Field::TRANSIENT, "\tTransient:  ");
      }
      total_faces += count;
      ++i;
    }
    if (summary) {
      OUTPUT << " Number of face sets                  =" << std::setw(9) << nss.size() << "\t";
      OUTPUT << " Length of face list        =" << std::setw(9) << total_faces << "\n";
    }
  }

  void info_elementsets(Ioss::Region &region, bool summary)
  {
    Ioss::ElementSetContainer      nss = region.get_elementsets();
    Ioss::ElementSetContainer::const_iterator i = nss.begin();
    int total_elements = 0;
    while (i != nss.end()) {
      int    count     = (*i)->get_property("entity_count").get_int();
      if (!summary) {
	OUTPUT << '\n' << (*i)->type_string() << " " << std::setw(16)  << (*i)->name()
	       << " id: " << std::setw(6) << id(*i)   << ", "
	       << std::setw(8) << count << " elements" << "\n";
	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
	info_fields(*i, Ioss::Field::TRANSIENT, "\tTransient:  ");
      }
      total_elements += count;
      ++i;
    }
    if (summary) {
      OUTPUT << " Number of element sets               =" << std::setw(9) << nss.size() << "\t";
      OUTPUT << " Length of element list     =" << std::setw(9) << total_elements << "\n";
    }
  }

  void info_commsets(Ioss::Region &region, bool summary)
  {
    Ioss::CommSetContainer      css = region.get_commsets();
    Ioss::CommSetContainer::const_iterator i = css.begin();
    while (i != css.end()) {
      std::string type      = (*i)->get_property("entity_type").get_string();
      ++i;
    }
    OUTPUT << '\n';
  }

  void info_fields(Ioss::GroupingEntity *ige,
		   Ioss::Field::RoleType role,
		   const std::string &header)
  {
    Ioss::NameList fields;
    ige->field_describe(role, &fields);

    if (fields.empty())
      return;
    
    if (!header.empty()) {
      OUTPUT << header;
    }
    // Iterate through results fields and transfer to output
    // database...  
    Ioss::NameList::const_iterator IF;
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;

      const Ioss::VariableType *var_type = ige->get_field(field_name).raw_storage();
      int comp_count = var_type->component_count();
      OUTPUT << std::setw(16) << field_name << ":" << comp_count << " ";
    }
    if (!header.empty()) {
      OUTPUT << "\n";
    }
  }

  void info_field_data(Ioss::GroupingEntity *ige,
		       Ioss::Field::RoleType role,
		       bool do_connectivity)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields;
    Ioss::NameList::const_iterator IF;
    ige->field_describe(role, &state_fields);

    for (IF = state_fields.begin(); IF != state_fields.end(); ++IF) {
      std::string field_name = *IF;
      // All of the 'Ioss::EntityBlock' derived classes have a
      // 'connectivity' field, but it is only interesting on the
      // Ioss::ElementBlock class. On the other classes, it just
      // generates overhead...
      if (!do_connectivity && field_name == "connectivity")
	continue;

      info_field_data_internal(ige, field_name);
    }
  }

  void info_field_data_internal(Ioss::GroupingEntity *ige, const std::string &field_name)
  {
#if 0
    size_t isize = ige->get_field(field_name).get_size();
    ige->get_field_data(field_name, &data[0], isize);
    OUTPUT << "Field: " << field_name << " has a data storage size of " << isize << " bytes.\n";
#endif
  }

  void info_properties(Ioss::GroupingEntity *ige)
  {
#if 0
    Ioss::NameList names;
    ige->property_describe(&names);

    // Iterate through properties and transfer to output database...
    Ioss::NameList::const_iterator I;
    for (I = names.begin(); I != names.end(); ++I) {
      OUTPUT << *I << ", ";
    }
#endif
  }

  void show_step(int istep, double time)
  {
    OUTPUT.setf(std::ios::scientific);
    OUTPUT.setf(std::ios::showpoint);
    OUTPUT << "     Time step " << std::setw(5) << istep
	   << " at time " << std::setprecision(5) << time << '\n';
  }
}
