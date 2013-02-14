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
#include <Ioss_SurfaceSplit.h>
#include <Ioss_Utils.h>
#include <init/Ionit_Initializer.h>
#include <stddef.h>
#include <stdlib.h>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "Ioss_CommSet.h"
#include "Ioss_DBUsage.h"
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
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_VariableType.h"

#include "info_interface.h"

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

  void show_step(int istep, double time);

  void info_nodeblock(Ioss::Region &region, const Info::Interface &interface, bool summary);
  void info_edgeblock(Ioss::Region &region, bool summary);
  void info_faceblock(Ioss::Region &region, bool summary);
  void info_elementblock(Ioss::Region &region, const Info::Interface &interface, bool summary);

  void info_nodesets(Ioss::Region &region, bool summary);
  void info_edgesets(Ioss::Region &region, bool summary);
  void info_facesets(Ioss::Region &region, bool summary);
  void info_elementsets(Ioss::Region &region, bool summary);

  void info_sidesets(Ioss::Region &region, const Info::Interface &interface, bool summary);
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
		 Info::Interface& interface);

  std::string name(Ioss::GroupingEntity *entity) {
    return entity->type_string() + " '" + entity->name() + "'";
  }

  int64_t id(Ioss::GroupingEntity *entity) {
    int64_t id = -1;
    if (entity->property_exists("id")) {
      id = entity->get_property("id").get_int();
    }
    return id;
  }


}
void hex_volume(Ioss::ElementBlock *block, const std::vector<double> &coordinates);

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
  
  Info::Interface interface;
  interface.parse_options(argc, argv);
  
  std::string in_type = "exodusII";

  codename = argv[0];
  size_t ind = codename.find_last_of("/", codename.size());
  if (ind != std::string::npos)
    codename = codename.substr(ind+1, codename.size());

  Ioss::Init::Initializer io;
#ifndef NO_XDMF_SUPPORT
  Ioxf::Initializer ioxf;
#endif

  OUTPUT << "Input:    '" << interface.filename()  << "', Type: " << interface.type()  << '\n';
  OUTPUT << '\n';

  file_info(interface.filename(), interface.type(), interface);

  OUTPUT << "\n" << codename << " execution successful.\n";
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return EXIT_SUCCESS;
}

namespace {
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

  void file_info(const std::string& inpfile, const std::string& input_type, Info::Interface& interface)
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

    dbi->set_surface_split_type(Ioss::int_to_surface_split(interface.surface_split_scheme()));
    dbi->set_field_separator(interface.field_suffix_separator());
    if (interface.ints_64_bit())
      dbi->set_int_byte_size_api(Ioss::USE_INT64_API);
    
    // NOTE: 'region' owns 'db' pointer at this time...
    Ioss::Region region(dbi, "region_1");

    // Get all properties of input database...
    bool summary = true;
    info_properties(&region);
    info_nodeblock(region,    interface, summary);
    info_edgeblock(region,    summary);
    info_faceblock(region,    summary);
    info_elementblock(region, interface, summary);

    info_nodesets(region,     summary);
    info_edgesets(region,     summary);
    info_facesets(region,     summary);
    info_elementsets(region,  summary);

    info_sidesets(region,     interface, summary);
    info_commsets(region,     summary);

    if (region.property_exists("state_count") && region.get_property("state_count").get_int() > 0) {
      std::pair<int, double> state_time_max = region.get_max_time();
      std::pair<int, double> state_time_min = region.get_min_time();
      OUTPUT << " Number of time steps on database     =" << std::setw(12)
	     << region.get_property("state_count").get_int() << "\n"
	     << "    Minimum time = " << state_time_min.second << " at step " << state_time_min.first << "\n"
	     << "    Maximum time = " << state_time_max.second << " at step " << state_time_max.first << "\n\n";
    }

    if (!interface.summary()) {
      summary = false;
      info_properties(&region);
      info_nodeblock(region,    interface, summary);
      info_edgeblock(region,    summary);
      info_faceblock(region,    summary);
      info_elementblock(region, interface, summary);
      
      info_nodesets(region,     summary);
      info_edgesets(region,     summary);
      info_facesets(region,     summary);
      info_elementsets(region,  summary);
      
      info_sidesets(region,     interface, summary);
      info_commsets(region,     summary);
    }
    
    if (interface.compute_volume()) {
      element_volume(region);
    }
  }


  void info_nodeblock(Ioss::Region &region, const Info::Interface &interface, bool summary)
  {
    Ioss::NodeBlockContainer    nbs = region.get_node_blocks();
    Ioss::NodeBlockContainer::const_iterator i = nbs.begin();
    while (i != nbs.end()) {
      //      std::string name      = (*i)->name();
      int64_t    num_nodes = (*i)->get_property("entity_count").get_int();
      int64_t    degree    = (*i)->get_property("component_degree").get_int();
      int64_t    num_attrib= (*i)->get_property("attribute_count").get_int();
      if (summary) {
	OUTPUT << " Number of spatial dimensions =" << std::setw(12) << degree << "\n";
	OUTPUT << " Number of nodeblocks         =" << std::setw(12) << 1 << "\t";
	OUTPUT << " Number of nodes            =" << std::setw(12) << num_nodes << "\n";
      } else {
	OUTPUT << '\n' << name(*i) 
	       << std::setw(12) << num_nodes << " nodes, "
	       << std::setw(3) << num_attrib << " attributes.\n";
	if (interface.check_node_status()) {
	  std::vector<char> node_status;
	  std::vector<int64_t>  ids;
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

  void info_elementblock(Ioss::Region &region, const Info::Interface &interface, bool summary)
  {
    Ioss::ElementBlockContainer ebs = region.get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
    int64_t total_elements = 0;
    while (i != ebs.end()) {
      int64_t    num_elem  = (*i)->get_property("entity_count").get_int();
      total_elements += num_elem;

      if (!summary) {
	std::string type      = (*i)->get_property("topology_type").get_string();
	int64_t    num_attrib= (*i)->get_property("attribute_count").get_int();
	OUTPUT << '\n' << name(*i)
	       << " id: " << std::setw(6) << id(*i)
	       << ", topology: " << std::setw(10) << type << ", "
	       << std::setw(12) << num_elem << " elements, "
	       << std::setw(3) << num_attrib << " attributes.";

	info_fields(*i, Ioss::Field::ATTRIBUTE, "\n\tAttributes: ");

	if (interface.adjacencies()) {
	  std::vector<std::string> blocks;
	  (*i)->get_block_adjacencies(blocks);
	  OUTPUT << "\n\tAdjacent to  " << blocks.size() << " element block(s):\t";
	  std::vector<std::string>::iterator b = blocks.begin();
	  while (b != blocks.end()) {
	    OUTPUT << *b++ << "  ";
	  }
	}
	info_fields(*i, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
	OUTPUT << "\n";
      }

      ++i;
    }
    if (summary) {
      OUTPUT << " Number of element blocks     =" << std::setw(12) << ebs.size() << "\t";
      OUTPUT << " Number of elements         =" << std::setw(12) << total_elements << "\n";
    }

  }

  void info_edgeblock(Ioss::Region &region, bool summary)
  {
    Ioss::EdgeBlockContainer ebs = region.get_edge_blocks();
    Ioss::EdgeBlockContainer::const_iterator i = ebs.begin();
    int64_t total_edges = 0;
    while (i != ebs.end()) {
      int64_t    num_edge  = (*i)->get_property("entity_count").get_int();
      total_edges += num_edge;

      if (!summary) {
	std::string type      = (*i)->get_property("topology_type").get_string();
	int64_t    num_attrib= (*i)->get_property("attribute_count").get_int();
	OUTPUT << '\n' << name(*i)
	       << " id: " << std::setw(6) << id(*i)
	       << ", topology: " << std::setw(10) << type << ", "
	       << std::setw(12) << num_edge << " edges, "
	       << std::setw(3) << num_attrib << " attributes.\n";

	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

#if 0
	std::vector<std::string> blocks;
	(*i)->get_block_adjacencies(blocks);
	OUTPUT << "\tAdjacent to  " << blocks.size() << " edge block(s):\t";
	std::vector<std::string>::iterator b = blocks.begin();
	while (b != blocks.end()) {
	  OUTPUT << *b++ << "  ";
	}
#endif
	info_fields(*i, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
	OUTPUT << "\n";
      }

      ++i;
    }
    if (summary) {
      OUTPUT << " Number of edge blocks        =" << std::setw(12) << ebs.size() << "\t";
      OUTPUT << " Number of edges            =" << std::setw(12) << total_edges << "\n";
    }

  }

  void info_faceblock(Ioss::Region &region, bool summary)
  {
    Ioss::FaceBlockContainer ebs = region.get_face_blocks();
    Ioss::FaceBlockContainer::const_iterator i = ebs.begin();
    int64_t total_faces = 0;
    while (i != ebs.end()) {
      int64_t    num_face  = (*i)->get_property("entity_count").get_int();
      total_faces += num_face;

      if (!summary) {
	std::string type      = (*i)->get_property("topology_type").get_string();
	int64_t    num_attrib= (*i)->get_property("attribute_count").get_int();
	OUTPUT << '\n' << name(*i)
	       << " id: " << std::setw(6) << id(*i)
	       << ", topology: " << std::setw(10) << type << ", "
	       << std::setw(12) << num_face << " faces, "
	       << std::setw(3) << num_attrib << " attributes.\n";

	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");

#if 0
	std::vector<std::string> blocks;
	(*i)->get_block_adjacencies(blocks);
	OUTPUT << "\tAdjacent to  " << blocks.size() << " face block(s):\t";
	std::vector<std::string>::iterator b = blocks.begin();
	while (b != blocks.end()) {
	  OUTPUT << *b++ << "  ";
	}
#endif
	info_fields(*i, Ioss::Field::TRANSIENT, "\n\tTransient:  ");
	OUTPUT << "\n";
      }

      ++i;
    }
    if (summary) {
      OUTPUT << " Number of face blocks        =" << std::setw(12) << ebs.size() << "\t";
      OUTPUT << " Number of faces            =" << std::setw(12) << total_faces << "\n";
    }

  }

  void info_sidesets(Ioss::Region &region, const Info::Interface &interface, bool summary)
  {
    Ioss::SideSetContainer      fss = region.get_sidesets();
    Ioss::SideSetContainer::const_iterator i = fss.begin();
    int64_t total_sides = 0;
    while (i != fss.end()) {
      if (!summary) {
	OUTPUT << '\n' << (*i)->type_string() << " " << std::setw(16)  << (*i)->name()
	       << " id: " << std::setw(6)<< id(*i) << ":";
	if (interface.adjacencies()) {
	  std::vector<std::string> blocks;
	  (*i)->block_membership(blocks);
	  OUTPUT << "\n\tTouches " << blocks.size() << " element block(s):\t";
	  std::vector<std::string>::iterator b = blocks.begin();
	  while (b != blocks.end()) {
	    OUTPUT << *b++ << "  ";
	  }
	  OUTPUT << "\n";
	}
      }
      if (!summary) {
	OUTPUT << "\n\tContains: \n";
      }
      
      Ioss::SideBlockContainer fbs = (*i)->get_side_blocks();
      Ioss::SideBlockContainer::const_iterator j = fbs.begin();
      while (j != fbs.end()) {
	int64_t    num_side  = (*j)->get_property("entity_count").get_int();
	if (!summary) {
	  std::string fbtype    = (*j)->get_property("topology_type").get_string();
	  std::string partype   = (*j)->get_property("parent_topology_type").get_string();
	  OUTPUT << "\t\t"
		 << (*j)->type_string() << " " << (*j)->name()
		 << ", "
		 << num_side << " " << fbtype << " sides"
		 << ", parent topology: " << partype 
		 << "\n";
	  if (interface.adjacencies()) {
	    std::vector<std::string> blocks;
	    (*j)->block_membership(blocks);
	    OUTPUT << "\t\t\tTouches " << blocks.size() << " element block(s):\t";
	    std::vector<std::string>::iterator b = blocks.begin();
	    while (b != blocks.end()) {
	      OUTPUT << *b++ << "  ";
	    }
	    OUTPUT << "\n";
	  }
	  OUTPUT << "\n";
	  info_fields(*j, Ioss::Field::ATTRIBUTE, "\t\tAttributes: ");
	  info_fields(*j, Ioss::Field::TRANSIENT, "\t\tTransient:  ");
	}
	total_sides += num_side;
	++j;
      }
      ++i;
    }

    if (summary) {
      OUTPUT << " Number of side sets          =" << std::setw(12) << fss.size() << "\t";
      OUTPUT << " Number of element sides    =" << std::setw(12) << total_sides << "\n";
    }
  }
  
  void info_nodesets(Ioss::Region &region, bool summary)
  {
    Ioss::NodeSetContainer      nss = region.get_nodesets();
    Ioss::NodeSetContainer::const_iterator i = nss.begin();
    int64_t total_nodes = 0;
    while (i != nss.end()) {
      int64_t    count     = (*i)->get_property("entity_count").get_int();
      int64_t    num_attrib= (*i)->get_property("attribute_count").get_int();
      int64_t    num_dist  = (*i)->get_property("distribution_factor_count").get_int();
      if (!summary) {
	OUTPUT << '\n' << (*i)->type_string() << " " << std::setw(16)  << (*i)->name()
	       << " id: " << std::setw(6) << id(*i)   << ", "
	       << std::setw(8) << count << " nodes" 
	       << std::setw(3) << num_attrib << " attributes"
	       << std::setw(8) << num_dist << " distribution factors.\n";
	info_fields(*i, Ioss::Field::ATTRIBUTE, "\tAttributes: ");
	info_fields(*i, Ioss::Field::TRANSIENT, "\tTransient:  ");
      }
      total_nodes += count;
      ++i;
    }
    if (summary) {
      OUTPUT << " Number of nodal point sets   =" << std::setw(12) << nss.size() << "\t";
      OUTPUT << " Length of node list        =" << std::setw(12) << total_nodes << "\n";
    }
  }

  void info_edgesets(Ioss::Region &region, bool summary)
  {
    Ioss::EdgeSetContainer      nss = region.get_edgesets();
    Ioss::EdgeSetContainer::const_iterator i = nss.begin();
    int64_t total_edges = 0;
    while (i != nss.end()) {
      int64_t    count     = (*i)->get_property("entity_count").get_int();
      int64_t    num_attrib= (*i)->get_property("attribute_count").get_int();
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
      OUTPUT << " Number of edge sets          =" << std::setw(12) << nss.size() << "\t";
      OUTPUT << " Length of edge list        =" << std::setw(12) << total_edges << "\n";
    }
  }

  void info_facesets(Ioss::Region &region, bool summary)
  {
    Ioss::FaceSetContainer      nss = region.get_facesets();
    Ioss::FaceSetContainer::const_iterator i = nss.begin();
    int64_t total_faces = 0;
    while (i != nss.end()) {
      int64_t    count     = (*i)->get_property("entity_count").get_int();
      int64_t    num_attrib= (*i)->get_property("attribute_count").get_int();
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
      OUTPUT << " Number of face sets          =" << std::setw(12) << nss.size() << "\t";
      OUTPUT << " Length of face list        =" << std::setw(12) << total_faces << "\n";
    }
  }

  void info_elementsets(Ioss::Region &region, bool summary)
  {
    Ioss::ElementSetContainer      nss = region.get_elementsets();
    Ioss::ElementSetContainer::const_iterator i = nss.begin();
    int64_t total_elements = 0;
    while (i != nss.end()) {
      int64_t    count     = (*i)->get_property("entity_count").get_int();
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
      OUTPUT << " Number of element sets       =" << std::setw(12) << nss.size() << "\t";
      OUTPUT << " Length of element list     =" << std::setw(12) << total_elements << "\n";
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
