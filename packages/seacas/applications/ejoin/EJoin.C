// Copyright(C) 2010 Sandia Corporation.
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
#include <iostream>
#include <iomanip>
#include <iterator>
#include <set>
#include <vector>
#include <numeric>
#include <map>
#include <string>
#include <exception>

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <cfloat>
#include <cmath>
#include <limits>

#include <cstring>
#include <ctime>
#include <sys/times.h>
#include <ctype.h>
#include <sys/utsname.h>

#include "adler.h"
#include "smart_assert.h"
#include "to_string.h"

#include <Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Transform.h>

#include "CodeTypes.h"
#include "SystemInterface.h"
#include "match_xyz.h"
#include "mapping.h"
#include "Vector3.h"
#include "Version.h"

namespace {
  bool valid_variable(const std::string variable, int id, const StringIdVector &variable_list);
  void define_global_fields(Ioss::Region &output_region, RegionVector &part_mesh,
			    const StringIdVector &variable_list);
  void define_nodal_fields(Ioss::Region &output_region, RegionVector &part_mesh,
			   const StringIdVector &variable_list,
			   SystemInterface &interface);
  void define_element_fields(Ioss::Region &output_region, RegionVector &part_mesh,
			     const StringIdVector &variable_list);
  void define_nset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
			  const StringIdVector &variable_list);
  void define_sset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
			  const StringIdVector &variable_list);
  void define_nodal_nodeset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
				   const StringIdVector &variable_list,
				   SystemInterface &interface);

  void output_nodeblock(Ioss::Region &output_region, RegionVector &part_mesh,
			const Map &local_node_map, Map &global_node_map);
  void output_elementblock(Ioss::Region &output_region, RegionVector &part_mesh,
			   const Map &local_node_map, const Map &local_element_map);
  void output_nodeset(Ioss::Region &output_region, RegionVector &part_mesh,
		      const Map &local_node_map);
  void output_sideset(Ioss::Region &output_region, RegionVector &part_mesh,
		      const Map& local_element_map);
  void output_nodal_nodeset(Ioss::Region &output_region, RegionVector &part_mesh,
			    SystemInterface &interface, const Map &local_node_map);
  void output_transient_state(Ioss::Region &output_region, RegionVector &part_mesh,
			      double time, const Map &local_node_map, SystemInterface &interface);
  void process_nset_omissions(RegionVector &part_mesh, const Omissions &omit);
  void process_sset_omissions(RegionVector &part_mesh, const Omissions &omit);

  template <typename T>
  bool approx_equal(T v1, T v2)
  {
#if 1
    static const T tolerance = 100.0 * std::numeric_limits<float>::epsilon();
    double d1 = std::fabs(v1 - v2);
    double d2 = std::fabs(v1+v2)*tolerance;
    return d1 <= d2;
#else
    return (float)v1 == (float)v2;
#endif
  }
}

extern void add_to_log(const char *name, int elapsed);
extern double ejoin_timer();

#include <exodusII.h>

namespace {
  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void create_nodal_nodeset(Ioss::Region &region, Ioss::Region &output_region, bool debug);
  void transfer_fields(Ioss::GroupingEntity *ige,
                       Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role,
                       const std::string &prefix = "");

  void transfer_field_data(Ioss::GroupingEntity *ige,
                           Ioss::GroupingEntity *oge,
                           Ioss::Field::RoleType role,
                           const std::string &prefix = "",
                           bool transfer_connectivity = true);

  void transfer_field_data_internal(Ioss::GroupingEntity *ige,
                                    Ioss::GroupingEntity *oge,
                                    const std::string &field_name);


  std::string time_stamp(const std::string &format);
}

std::string tsFormat = "[%H:%M:%S] ";

// prototypes

int ejoin(SystemInterface &interface, std::vector<Ioss::Region*> &part_mesh);

unsigned int debug_level = 0;

int main(int argc, char* argv[])
{
#if defined(__LIBCATAMOUNT__)
  setlinebuf(stderr);
#endif
  try {
  time_t begin_time = time(NULL);
  SystemInterface::show_version();
  Ioss::Init::Initializer io;

  SystemInterface interface;
  bool ok = interface.parse_options(argc, argv);

  if (!ok) {
    std::cerr << "\nERROR: Problems parsing command line arguments.\n\n";
    exit(EXIT_FAILURE);
  }

  debug_level = interface.debug();

  if (debug_level & 64)
    ex_opts(EX_VERBOSE|EX_DEBUG);
  else
    ex_opts(0);

  int error = 0;

  const Omissions &omissions = interface.block_omissions();
  std::vector<Ioss::Region*> part_mesh(interface.inputFiles_.size());
  for (size_t p=0; p < interface.inputFiles_.size(); p++) {
    Ioss::DatabaseIO *dbi = Ioss::IOFactory::create("exodusII", interface.inputFiles_[p],
						    Ioss::READ_RESTART, (MPI_Comm)MPI_COMM_WORLD);
    if (dbi == NULL || !dbi->ok(true))
      std::exit(EXIT_FAILURE);

    dbi->set_surface_split_type(Ioss::SPLIT_BY_DONT_SPLIT);

    if (interface.disable_field_recognition())
      dbi->set_field_separator(0);
    
    if (!omissions[p].empty()) {
      dbi->set_block_omissions(omissions[p]);
    }

    // Generate a name for the region based on the part number...
    std::string name = "p" + to_string(p+1);
    // NOTE: region owns database pointer at this time...
    part_mesh[p] = new Ioss::Region(dbi, name);
    part_mesh[p]->property_add(Ioss::Property("block_omission_count", (int)omissions[p].size()));

    Vector3 offset = interface.offset();
    if (p > 0 && (offset.x != 0.0 || offset.y != 0.0 || offset.z != 0.0)) {
      Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
      Ioss::Field coord = nb->get_field("mesh_model_coordinates");
      Ioss::Transform* transform = Iotr::Factory::create("offset3D");
      assert(transform != NULL);
      std::vector<double> values(3);
      values[0] = offset.x * p; values[1] = offset.y * p; values[2] = offset.z * p;
      transform->set_properties("offset", values);
      coord.add_transform(transform);
      nb->field_erase("mesh_model_coordinates");
      nb->field_add(coord);
    }
  }

  process_nset_omissions(part_mesh, interface.nset_omissions());
  process_sset_omissions(part_mesh, interface.sset_omissions());

  ejoin(interface, part_mesh);
  
  for (size_t p=0; p < interface.inputFiles_.size(); p++) {
    delete part_mesh[p];
  }
  
  time_t end_time = time(NULL);
  add_to_log(argv[0], (int)(end_time-begin_time));
  return (error);
  }
  catch (std::exception &e) {
    std::cerr << "ERROR: Standard exception: " << e.what() << std::endl;
  }
}

int ejoin(SystemInterface &interface, std::vector<Ioss::Region*> &part_mesh)
{
  size_t part_count = interface.inputFiles_.size();
  SMART_ASSERT(part_count == part_mesh.size());
  
  Ioss::DatabaseIO *dbo = Ioss::IOFactory::create("exodusII", interface.outputName_, Ioss::WRITE_RESTART,
						  (MPI_Comm)MPI_COMM_WORLD);
  if (dbo == NULL || !dbo->ok(true)) {
    std::exit(EXIT_FAILURE);
  }
  
  // NOTE: 'output_region' owns 'dbo' pointer at this time
  Ioss::Region output_region(dbo, "ejoin_output_region");

  output_region.begin_mode(Ioss::STATE_DEFINE_MODEL);
  
  output_region.property_add(Ioss::Property("code_name", qainfo[0]));
  output_region.property_add(Ioss::Property("code_version", qainfo[2]));

  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);

  int node_offset = 0;
  int element_offset = 0;
  for (size_t p = 0; p < part_count; p++) {
    part_mesh[p]->property_add(Ioss::Property("node_offset", node_offset));
    part_mesh[p]->property_add(Ioss::Property("element_offset", element_offset));
    int local_node_count = part_mesh[p]->get_property("node_count").get_int();
    int local_elem_count = part_mesh[p]->get_property("element_count").get_int();
    node_offset += local_node_count;
    element_offset += local_elem_count;
  }

  int node_count = node_offset; // Sum of nodes in part meshes.
  Map local_node_map(node_count);
  Map global_node_map;

  // This is the map from local element position to global element
  // position (0-based). If there are no element block omissions, then
  // the map is simply [0..number_elements). If there are omissions,
  // Then local_element_map[j] will be -1 for an omitted element.
  // If not omitted, local_element_map[part_offset+j] gives the global
  // position of local element j in the current part.
  Map local_element_map(element_offset);
  build_local_element_map(part_mesh, local_element_map);
			  
  // Need a map from local node to global node.  Position in the map
  // is part_mesh.offset+local_node_position Return value is position
  // in the global node list.  If no mapping, then the list is simply
  // 0..number_of_global_nodes where number_of_global_nodes is the sum
  // of the node counts in the individual files.
  if (interface.match_node_ids()) {
    build_reverse_node_map(output_region, part_mesh, global_node_map, local_node_map);
  } else if (interface.match_node_xyz()) {
    match_node_xyz(part_mesh, interface.tolerance(), global_node_map, local_node_map);
  } else if (!interface.block_omissions().empty()) {
    // At least 1 element block has been omitted in at least 1 part.
    // Eliminate all nodes that were only connected to the omitted element blocks.
    eliminate_omitted_nodes(part_mesh, global_node_map, local_node_map);
  } else {
    global_node_map.resize(local_node_map.size());
    for (size_t i=0; i < local_node_map.size(); i++) {
      local_node_map[i] = i;
      global_node_map[i] = i+1;
    }
  }
	     
  node_count = global_node_map.size();
  size_t merged = local_node_map.size() - global_node_map.size();
  if (merged > 0) {
    std::cerr << "*** " << merged << " Nodes were merged/omitted.\n";
  }

  // Verify nodemap...
  std::vector<int> glob(node_count);
  for (size_t i=0; i<local_node_map.size(); i++) {
    if (local_node_map[i] >= 0)
      glob[local_node_map[i]] = 1;
  }
  for (size_t i=0; i<glob.size(); i++) {
    SMART_ASSERT(glob[i] == 1);
  }
  // Transfer some common data...
  output_region.property_add(part_mesh[0]->get_property("title"));
  output_region.property_add(part_mesh[0]->get_property("spatial_dimension"));
  int spatial_dimension = part_mesh[0]->get_property("spatial_dimension").get_int();

  // Define a node block...  
  std::string block_name = "nodeblock_1";
  Ioss::NodeBlock *block = new Ioss::NodeBlock(output_region.get_database(), block_name,
					       node_count, spatial_dimension);
  block->property_add(Ioss::Property("id", 1));
  
  output_region.add(block);

  // Add element blocks, nodesets, sidesets
  for (size_t p = 0; p < part_count; p++) {
    transfer_elementblock(*part_mesh[p], output_region, false);
    if (interface.convert_nodes_to_nodesets(p+1)) {
      create_nodal_nodeset(*part_mesh[p], output_region, false);
    }
    if (!interface.omit_nodesets()) {
      transfer_nodesets(*part_mesh[p], output_region, false);
    }
    if (!interface.omit_sidesets()) {
      transfer_sidesets(*part_mesh[p], output_region, false);
    }
  }
    
  output_region.end_mode(Ioss::STATE_DEFINE_MODEL);

  output_region.begin_mode(Ioss::STATE_MODEL);

  output_nodeblock(output_region, part_mesh, local_node_map, global_node_map);
  output_elementblock(output_region, part_mesh, local_node_map, local_element_map);
  output_nodal_nodeset(output_region, part_mesh, interface, local_node_map);

  if (!interface.omit_nodesets())
    output_nodeset(output_region, part_mesh, local_node_map);
  if (!interface.omit_sidesets()) {
    output_sideset(output_region, part_mesh, local_element_map);
  }

  output_region.end_mode(Ioss::STATE_MODEL);
  
  // ####################TRANSIENT DATA SECTION###########################
  // ***********************************************************************
  // 9. Get Variable Information and names

  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);

  output_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);
  
  define_global_fields(output_region, part_mesh, interface.global_var_names());

  define_nodal_fields(output_region, part_mesh, interface.node_var_names(), interface);
  define_nodal_nodeset_fields(output_region, part_mesh, interface.node_var_names(), interface);

  define_element_fields(output_region, part_mesh, interface.elem_var_names());
  
  if (!interface.omit_nodesets())
    define_nset_fields(output_region, part_mesh, interface.nset_var_names());
  if (!interface.omit_sidesets()) {
    define_sset_fields(output_region, part_mesh, interface.sset_var_names());
  }

  output_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);

  // Get database times...
  // Different parts can have either no times or the times must match
  //! \todo If ts_min, ts_max, ts_step is specified, then only check steps in that range...
  std::vector<double> global_times;
  
  for (size_t p = 0; p < part_count; p++) {
    int nts = part_mesh[p]->get_property("state_count").get_int();
    std::vector<double> times(nts);
    // If multiple databases have timesteps, they must all match.
    for (int i=0; i<nts; i++) {
      times[i] = part_mesh[p]->get_state_time(i+1);
    }

    if (nts > 0) {
      if (global_times.empty()) {
	std::copy(times.begin(), times.end(), std::back_inserter(global_times));
      } else {
	if (global_times.size() != times.size()) {
	  std::cerr << "Time step sizes must match.";
	  SMART_ASSERT(global_times.size() == times.size())(global_times.size())(times.size())(p);
	  exit(EXIT_FAILURE);
	}

	for (size_t i=0; i < global_times.size(); i++) {
	  if (!approx_equal(global_times[i], times[i])) {
	    std::cerr << "Time step " << i << " in part " << p+1
		      << " does not match time steps in previous part(s): previous: "
		      << global_times[i] << ", current: " << times[i] << "\n";
	    exit(EXIT_FAILURE);
	  }
	}
      }
    }
  }

  output_region.begin_mode(Ioss::STATE_TRANSIENT);
  std::cerr << "\n";
  int ts_min  = interface.step_min();
  int ts_max  = interface.step_max();
  int ts_step = interface.step_interval();
  int num_steps = (int)global_times.size();

  if (ts_min == -1 && ts_max == -1) {
    ts_min = num_steps;
    ts_max = num_steps;
  }
  ts_max = ts_max < num_steps ? ts_max : num_steps;

  for (int step=ts_min-1; step < ts_max; step+=ts_step) {
    int ostep = output_region.add_state(global_times[step]);
    output_region.begin_state(ostep);
    output_transient_state(output_region, part_mesh, global_times[step], local_node_map, interface);
    std::cerr << "\rWrote step " << std::setw(4) << step+1 << "/" << global_times.size() << ", time "
	      << std::scientific << std::setprecision(4) << global_times[step];
    output_region.end_state(ostep);
  }
  std::cerr << "\n";
  output_region.end_mode(Ioss::STATE_TRANSIENT);
  
  /*************************************************************************/
  // EXIT program
  if (debug_level & 1)
    std::cerr << time_stamp(tsFormat);
  output_region.output_summary(std::cout);
  std::cerr << "******* END *******\n";
  return(0);
}

namespace {
  bool entity_is_omitted(Ioss::GroupingEntity *block) {
    bool omitted = false;
    if (block->property_exists("omitted"))
      omitted = (block->get_property("omitted").get_int() == 1);
    return omitted;
  }

  void transfer_elementblock(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    static int used_blocks = 0;
    std::string prefix = region.name();
    
    Ioss::ElementBlockContainer ebs = region.get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator i = ebs.begin();
    int total_elements = 0;
    while (i != ebs.end()) {
      Ioss::ElementBlock *eb = *i;
      if (!entity_is_omitted(eb)) {
	std::string name      = prefix + "_" + eb->name();
	if (output_region.get_element_block(name) != NULL) {
	  std::cerr << "ERROR: Duplicate element blocks named '" << name << "'\n";
	  exit(EXIT_FAILURE);
	}
	if (debug) std::cerr << name << ", ";
	std::string type      = eb->get_property("topology_type").get_string();
	int    num_elem  = eb->get_property("entity_count").get_int();
	total_elements += num_elem;

	if (num_elem > 0) {
	  Ioss::ElementBlock *ebn = new Ioss::ElementBlock(output_region.get_database(), name, type,
							   num_elem);
	  ebn->property_add(Ioss::Property("original_block_order", used_blocks++));
	  output_region.add(ebn);
	  transfer_fields(eb, ebn, Ioss::Field::ATTRIBUTE);

	  if (eb->property_exists("original_topology_type")) {
	    std::string oes = eb->get_property("original_topology_type").get_string();

	    // Set the new property
	    ebn->property_add(Ioss::Property("original_topology_type", oes));
	  }
	}
      }
      ++i;
    }
  }

  void transfer_sidesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    std::string prefix = region.name();

    Ioss::SideSetContainer      fss = region.get_sidesets();
    Ioss::SideSetContainer::const_iterator i = fss.begin();
    int total_sides = 0;
    while (i != fss.end()) {
      if (!entity_is_omitted(*i)) {
	std::string name      = prefix + "_" + (*i)->name();
	if (output_region.get_sideset(name) != NULL) {
	  std::cerr << "ERROR: Duplicate side sets named '" << name << "'\n";
	  exit(EXIT_FAILURE);
	}
	if (debug) std::cerr << name << ", ";
	Ioss::SideSet *surf = new Ioss::SideSet(output_region.get_database(), name);

	Ioss::SideBlockContainer fbs = (*i)->get_side_blocks();
	Ioss::SideBlockContainer::const_iterator j = fbs.begin();
	while (j != fbs.end()) {
	  std::string fbname    = prefix + "_" + (*j)->name();
	  if (debug) std::cerr << fbname << ", ";
	  std::string fbtype    = (*j)->get_property("topology_type").get_string();
	  std::string partype   = (*j)->get_property("parent_topology_type").get_string();
	  int    num_side  = (*j)->get_property("entity_count").get_int();
	  total_sides += num_side;

	  Ioss::SideBlock *block = new Ioss::SideBlock(output_region.get_database(), fbname, fbtype,
						       partype, num_side);
	  surf->add(block);
	  ++j;
	}
	output_region.add(surf);
      }
      ++i;
    }
  }

  // Create a nodeset on the output region consisting of all the nodes
  // in the input region.
  void create_nodal_nodeset(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    std::string prefix = region.name();

    std::string name = prefix + "_nodes";
    if (output_region.get_nodeset(name) != NULL) {
      std::cerr << "ERROR: Duplicate node sets named '" << name << "'\n";
      exit(EXIT_FAILURE);
    }
    if (debug) std::cerr << name << ", ";

    int count = region.get_property("node_count").get_int();
    Ioss::NodeSet *ns = new Ioss::NodeSet(output_region.get_database(), name, count);
    output_region.add(ns);
  }

  // Output the bulk data for a nodeset on the output region
  // consisting of all the nodes in the input region.
  void output_nodal_nodeset(Ioss::Region &output_region, RegionVector &part_mesh, SystemInterface &interface,
			    const Map &local_node_map)
  {
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      if (interface.convert_nodes_to_nodesets(p+1)) {
	std::string name = part_mesh[p]->name() + "_nodes";
	Ioss::NodeSet *ons = output_region.get_nodeset(name);
	SMART_ASSERT(ons != NULL)(name);

	Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
	SMART_ASSERT(nb != NULL);
	int node_offset = part_mesh[p]->get_property("node_offset").get_int();

	std::vector<int> nodelist;
	nb->get_field_data("ids", nodelist);
	for (size_t i=0; i < nodelist.size(); i++) {
	  int loc_node = part_mesh[p]->node_global_to_local(nodelist[i], true)-1;
	  int gpos = local_node_map[node_offset+loc_node];
	  if (gpos >= 0)
	    nodelist[i] = gpos+1;
	}
	ons->put_field_data("ids", nodelist);
      }
    }
  }
  
  void define_nodal_nodeset_fields(Ioss::Region &output_region, RegionVector &part_mesh, const StringIdVector &variable_list,
				   SystemInterface &interface)
  {
    if (!variable_list.empty() && variable_list[0].first == "none")
      return;

    size_t part_count = part_mesh.size();
    for (size_t p=0; p<part_count; p++) {

      if (interface.convert_nodes_to_nodesets(p+1)) {
	// Find nodeset in output region corresponding to the nodes in this part...
	std::string name = part_mesh[p]->name() + "_nodes";
	Ioss::NodeSet *ons = output_region.get_nodeset(name);
	SMART_ASSERT(ons != NULL)(name);
	
	Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
	SMART_ASSERT(nb != NULL);

	SMART_ASSERT(part_mesh[p]->get_property("node_count").get_int() == nb->get_property("entity_count").get_int());

	Ioss::NameList fields;
	nb->field_describe(Ioss::Field::TRANSIENT, &fields);
	Ioss::NameList::const_iterator IF;
	for (IF = fields.begin(); IF != fields.end(); ++IF) {
	  if (valid_variable(*IF, 0, variable_list)) {
	    Ioss::Field field = nb->get_field(*IF);
	    ons->field_add(field);
	  }
	}
      }
    }
  }

  void transfer_nodesets(Ioss::Region &region, Ioss::Region &output_region, bool debug)
  {
    std::string prefix = region.name();

    Ioss::NodeSetContainer      nss = region.get_nodesets();
    Ioss::NodeSetContainer::const_iterator i = nss.begin();
    while (i != nss.end()) {
      if (!entity_is_omitted(*i)) {
	std::string name = prefix + "_" + (*i)->name();
	if (output_region.get_nodeset(name) != NULL) {
	  std::cerr << "ERROR: Duplicate node sets named '" << name << "'\n";
	  exit(EXIT_FAILURE);
	}
	if (debug) std::cerr << name << ", ";
	int    count     = (*i)->get_property("entity_count").get_int();
	Ioss::NodeSet *ns = new Ioss::NodeSet(output_region.get_database(), name, count);
	output_region.add(ns);
	int id = (*i)->get_property("id").get_int();
	ns->property_add(Ioss::Property("id", id));
      }
      ++i;
    }
  }

  std::string time_stamp(const std::string &format)
  {
    if (format == "") {
      return std::string("");
    } else {
      const int length=256;
      static char time_string[length];

      time_t calendar_time = time(NULL);
      struct tm *local_time = localtime(&calendar_time);

      int error = strftime(time_string, length, format.c_str(), local_time);
      if (error != 0) {
        time_string[length-1] = (char)NULL;
        return std::string(time_string);
      } else {
        return std::string("[ERROR]");
      }
    }
  }

  template <typename T>
  void map_element_vars(int loffset, int goffset, int entity_count, 
			std::vector<T> &values,
			std::vector<T> &global_values,
			int *part_loc_elem_to_global)
  {
    // copy values to master element value information
    T* local_values = &values[0];
    for (int j = 0; j < entity_count; j++) {
      int global_block_pos = part_loc_elem_to_global[(j + loffset)] - goffset;
      global_values[global_block_pos] = local_values[j];
    }
  }

  template <typename T>
  void map_sideset_vars(int loffset, int entity_count, 
			std::vector<T> &values,
			std::vector<T> &global_values)
  {
    // copy values to master sideset value information
    T* local_values = &values[0];
    for (int j = 0; j < entity_count; j++) {
      global_values[j+loffset] = local_values[j];
    }
  }

  template <typename T, typename U>
  void map_nodeset_vars(U &, int, int, std::vector<T> &, std::vector<T> &)
  {
    SMART_ASSERT(1==0 && "Internal Error!");
  }
}

namespace {
  void output_nodeblock(Ioss::Region &output_region, RegionVector &part_mesh, const Map &local_node_map,
			Map &global_node_map)
  {
    Ioss::NodeBlock *onb = output_region.get_node_blocks()[0];
    SMART_ASSERT(onb != NULL);

    onb->put_field_data("ids", global_node_map);

    int spatial_dimension = output_region.get_property("spatial_dimension").get_int();
    std::vector<double> coord(global_node_map.size()*spatial_dimension);
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
      SMART_ASSERT(nb != NULL);
      std::vector<double> coordinates;
      nb->get_field_data("mesh_model_coordinates", coordinates);
      size_t node_count = nb->get_property("entity_count").get_int();
      int offset = part_mesh[p]->get_property("node_offset").get_int();
      for (size_t i=0; i < node_count; i++) {
	int glob_pos = local_node_map[i+offset];
	if (glob_pos >= 0) {
	  coord[glob_pos*spatial_dimension + 0] = coordinates[i*spatial_dimension + 0];
	  coord[glob_pos*spatial_dimension + 1] = coordinates[i*spatial_dimension + 1];
	  coord[glob_pos*spatial_dimension + 2] = coordinates[i*spatial_dimension + 2];
	}
      }
    }
    onb->put_field_data("mesh_model_coordinates", coord);
  }

  void output_elementblock(Ioss::Region &output_region, RegionVector &part_mesh,
			   const Map &local_node_map, const Map &local_element_map)
  {

    Ioss::ElementBlockContainer ebs = output_region.get_element_blocks();

    int element_count = output_region.get_property("element_count").get_int();
    Map ids(element_count);

    // Try to maintain the original element ids if possible...
    generate_element_ids(part_mesh, local_element_map, ids);

    int element_offset = 0;
    Ioss::ElementBlockContainer::const_iterator I = ebs.begin();
    while (I != ebs.end()) {
      (*I)->put_field_data("ids", &ids[element_offset], ids.size()*sizeof(int));
      element_offset += (*I)->get_property("entity_count").get_int();
      ++I;
    }

    SMART_ASSERT(element_offset == element_count);

    // Connectivity...
    I = ebs.begin();
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      Ioss::ElementBlockContainer iebs = part_mesh[p]->get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator J = iebs.begin();
      int node_offset = part_mesh[p]->get_property("node_offset").get_int();
      
      while (J != iebs.end()) {
	Ioss::ElementBlock *ieb = *J;
	std::string name = part_mesh[p]->name() + "_" + ieb->name();
	Ioss::ElementBlock *oeb = output_region.get_element_block(name);
	if (oeb != NULL) {
	  std::vector<int> connectivity;
	  ieb->get_field_data("connectivity", connectivity);
	
	  SMART_ASSERT(ieb->get_property("entity_count").get_int() == oeb->get_property("entity_count").get_int());
	  for (size_t i=0; i < connectivity.size(); i++) {

	    // connectivity is in part-global node ids.
	    // loc_node = the position of part-global node 'connectivity[i]' in the local [0..num_node)
	    // local_node_map[node_offset+loc_node] gives the position of this node in the global list
	    int loc_node = part_mesh[p]->node_global_to_local(connectivity[i], true)-1;
	    SMART_ASSERT(node_offset+loc_node < (int)local_node_map.size());
	    int gpos = local_node_map[node_offset+loc_node];
	    if (gpos >= 0) 
	      connectivity[i] = gpos+1;
	  }
	  oeb->put_field_data("connectivity", connectivity);
	  transfer_field_data(ieb, oeb, Ioss::Field::ATTRIBUTE);
	}
	++J;
      }
    }
  }

  void output_nodeset(Ioss::Region &output_region, RegionVector &part_mesh, const Map &local_node_map)
  {
    if (output_region.get_nodesets().empty())
      return;

    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      int node_offset = part_mesh[p]->get_property("node_offset").get_int();
      Ioss::NodeSetContainer ins = part_mesh[p]->get_nodesets();
      Ioss::NodeSetContainer::const_iterator J = ins.begin();
      while (J != ins.end()) {
	if (!entity_is_omitted(*J)) {
	  std::vector<int> nodelist;
	  (*J)->get_field_data("ids", nodelist);
	
	  std::string name = part_mesh[p]->name() + "_" + (*J)->name();
	  Ioss::NodeSet *ons = output_region.get_nodeset(name);
	  SMART_ASSERT(ons != NULL)(name);
	  SMART_ASSERT((*J)->get_property("entity_count").get_int() == ons->get_property("entity_count").get_int());

	  // This needs to make sure that the nodelist comes back as local id (1..numnodes)
	  for (size_t i=0; i < nodelist.size(); i++) {
	    int loc_node = part_mesh[p]->node_global_to_local(nodelist[i], true)-1;
	    int gpos = local_node_map[node_offset+loc_node];
	    if (gpos >= 0) 
	      nodelist[i] = gpos+1;
	  }
	  ons->put_field_data("ids", nodelist);
	}
	++J;
      }
    }
  }
  
  void output_sideset(Ioss::Region &output_region, RegionVector &part_mesh,
		      const Map& local_element_map)
  {
    Ioss::SideSetContainer os = output_region.get_sidesets();
    Ioss::SideSetContainer::const_iterator I = os.begin();

    Ioss::SideBlockContainer out_eb;
    // Put all output side blocks in the same list...
    while (I != os.end()) {
      Ioss::SideBlockContainer obs = (*I++)->get_side_blocks();
      std::copy(obs.begin(), obs.end(), std::back_inserter(out_eb));
    }
    
    // Assuming (with checks) that the output side blocks will be
    // iterated in same order as input side blocks...
    Ioss::SideBlockContainer::const_iterator II = out_eb.begin();

    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      int element_offset = part_mesh[p]->get_property("element_offset").get_int();

      Ioss::SideSetContainer is = part_mesh[p]->get_sidesets();
      Ioss::SideSetContainer::const_iterator J = is.begin();
      while (J != is.end()) {
	if (!entity_is_omitted(*J)) {
	  Ioss::SideBlockContainer ebs = (*J)->get_side_blocks();
	  Ioss::SideBlockContainer::const_iterator JJ = ebs.begin();

	  while (JJ != ebs.end()) {
	    SMART_ASSERT(part_mesh[p]->name() + "_" + (*JJ)->name() == (*II)->name())((*JJ)->name())((*II)->name());
	    SMART_ASSERT((*JJ)->get_property("entity_count").get_int() == (*II)->get_property("entity_count").get_int());
	    std::vector<int> elem_side_list;
	    (*JJ)->get_field_data("element_side_raw", elem_side_list);

	    // The 'elem_side_list' contains
	    // (local_element_position,side_ordinal) pairs. The
	    // 'local_element_position' is 1-based offset in the
	    // current part.  Need to map to its location in the
	    // output region...
	    for (size_t i=0; i < elem_side_list.size(); i+=2) { // just get the elem part of the pair...
	      int local_position = elem_side_list[i]-1;
	      int gpos = local_element_map[element_offset + local_position];
	      SMART_ASSERT(gpos >= 0)(gpos)(i); // Inactive elements should be filtered by Ioss
	      elem_side_list[i] = gpos+1;
	    }
	    (*II)->put_field_data("element_side_raw", elem_side_list);
	    ++JJ; ++II;
	  }
	}
	++J;
      }
    }
  }

  void output_globals(Ioss::Region &output_region, RegionVector &part_mesh,
		      double time, const IntVector &steps)
  {
    size_t part_count = part_mesh.size();
    for (size_t p=0; p<part_count; p++) {
      Ioss::NameList fields;
      part_mesh[p]->field_describe(Ioss::Field::TRANSIENT, &fields);
      Ioss::NameList::const_iterator IF;
      for (IF = fields.begin(); IF != fields.end(); ++IF) {
	std::vector<double> data;
	part_mesh[p]->get_field_data(*IF, data);
	output_region.put_field_data(*IF, data);
      }
    }
  }

  void output_nodal(Ioss::Region &output_region, RegionVector &part_mesh,
		    double time, const IntVector &steps, const Map &local_node_map,
		    SystemInterface &interface)
  {
    size_t part_count = part_mesh.size();

    Ioss::NodeBlock *onb = output_region.get_node_blocks()[0];
    SMART_ASSERT(onb != NULL);
    int node_count = onb->get_property("entity_count").get_int();

    Ioss::NameList fields;
    onb->field_describe(Ioss::Field::TRANSIENT, &fields);
    Ioss::NameList::const_iterator IF;
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      size_t comp_count = onb->get_field(*IF).raw_storage()->component_count();
      std::vector<double> data(node_count*comp_count);
      for (size_t p=0; p<part_count; p++) {
	if (!interface.convert_nodes_to_nodesets(p+1)) {
	  int offset = part_mesh[p]->get_property("node_offset").get_int();
	  Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
	  SMART_ASSERT(nb != NULL);
	  if (nb->field_exists(*IF)) {
	    SMART_ASSERT((int)comp_count == nb->get_field(*IF).raw_storage()->component_count());
	    std::vector<double> loc_data;
	    nb->get_field_data(*IF, loc_data);
	    size_t nc = nb->get_property("entity_count").get_int();
	    SMART_ASSERT(loc_data.size() == nc * comp_count);
	    for (size_t i=0; i < nc; i++) {
	      int glob_pos = local_node_map[offset+i];
	      if (glob_pos >= 0) {
		for (size_t j=0; j<comp_count; j++) {
		  data[glob_pos*comp_count+j] = loc_data[i*comp_count+j];
		}
	      }
	    }
	  }
	}
      }
      onb->put_field_data(*IF, data);
    }
  }

  void output_nodal_nodeset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
				   double time, const IntVector &steps, 
				   SystemInterface &interface)
  {
    size_t part_count = part_mesh.size();
    for (size_t p=0; p<part_count; p++) {
      if (interface.convert_nodes_to_nodesets(p+1)) {

	// Find nodeset in output region corresponding to the nodes in this part...
	std::string name = part_mesh[p]->name() + "_nodes";
	Ioss::NodeSet *ons = output_region.get_nodeset(name);
	SMART_ASSERT(ons != NULL)(name);
	
	Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
	SMART_ASSERT(nb != NULL);

	SMART_ASSERT(part_mesh[p]->get_property("node_count").get_int() == nb->get_property("entity_count").get_int());

	// NOTE: The node order in the output nodeset 'ons' was
	// defined as the same node order in the input nodeblock 'nb',
	// so we shouldn't have to do any reordering of the data at
	// this time--just read then write.
	Ioss::NameList fields;
	ons->field_describe(Ioss::Field::TRANSIENT, &fields);
	Ioss::NameList::const_iterator IF;
	std::vector<double> data;
	for (IF = fields.begin(); IF != fields.end(); ++IF) {
	  nb->get_field_data(*IF, data);
	  ons->put_field_data(*IF, data);
	}
      }
    }
  }

  void output_element(Ioss::Region &output_region, RegionVector &part_mesh,
		      double time, const IntVector &steps)
  {
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      Ioss::ElementBlockContainer iebs = part_mesh[p]->get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator J = iebs.begin();
      while (J != iebs.end()) {
	Ioss::ElementBlock *ieb = *J;
	std::string name = part_mesh[p]->name() + "_" + ieb->name();

	Ioss::ElementBlock *oeb = output_region.get_element_block(name);
	if (oeb != NULL) {
	  Ioss::NameList fields;
	  ieb->field_describe(Ioss::Field::TRANSIENT, &fields);
	  Ioss::NameList::const_iterator IF;
	  for (IF = fields.begin(); IF != fields.end(); ++IF) {
	    if (oeb->field_exists(*IF)) {
	      transfer_field_data_internal(ieb, oeb, *IF);
	    }
	  }
	}
	++J;
      }
    }
  }

  void output_nset(Ioss::Region &output_region, RegionVector &part_mesh,
		   double time, const IntVector &steps)
  {
    if (output_region.get_nodesets().empty())
      return;

    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      Ioss::NodeSetContainer ins = part_mesh[p]->get_nodesets();
      Ioss::NodeSetContainer::const_iterator J = ins.begin();
      while (J != ins.end()) {
	if (!entity_is_omitted(*J)) {
	  std::string name = part_mesh[p]->name() + "_" + (*J)->name();
	  Ioss::NodeSet *ons = output_region.get_nodeset(name);
	  SMART_ASSERT(ons != NULL)(name);

	  Ioss::NameList fields;
	  (*J)->field_describe(Ioss::Field::TRANSIENT, &fields);
	  Ioss::NameList::const_iterator IF;
	  for (IF = fields.begin(); IF != fields.end(); ++IF) {
	    if (ons->field_exists(*IF)) {
	      transfer_field_data_internal(*J, ons, *IF);
	    }
	  }
	}
	++J;
      }
    }
  }

  void output_sset(Ioss::Region &output_region, RegionVector &part_mesh,
		   double time, const IntVector &steps)
  {
    Ioss::SideSetContainer os = output_region.get_sidesets();
    if (os.empty())
      return;
    Ioss::SideSetContainer::const_iterator I = os.begin();

    Ioss::SideBlockContainer out_eb;
    // Put all output side blocks in the same list...
    while (I != os.end()) {
      Ioss::SideBlockContainer obs = (*I++)->get_side_blocks();
      std::copy(obs.begin(), obs.end(), std::back_inserter(out_eb));
    }
    
    // Assuming (with checks) that the output side blocks will be
    // iterated in same order as input side blocks...
    Ioss::SideBlockContainer::const_iterator II = out_eb.begin();

    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      Ioss::SideSetContainer is = part_mesh[p]->get_sidesets();
      Ioss::SideSetContainer::const_iterator J = is.begin();
      while (J != is.end()) {
	if (!entity_is_omitted(*J)) {
	  Ioss::SideBlockContainer ebs = (*J)->get_side_blocks();
	  Ioss::SideBlockContainer::const_iterator JJ = ebs.begin();

	  while (JJ != ebs.end()) {
	    SMART_ASSERT(part_mesh[p]->name() + "_" + (*JJ)->name() == (*II)->name());
	    Ioss::NameList fields;
	    (*JJ)->field_describe(Ioss::Field::TRANSIENT, &fields);
	    Ioss::NameList::const_iterator IF;
	    for (IF = fields.begin(); IF != fields.end(); ++IF) {
	      if ((*II)->field_exists(*IF)) {
		transfer_field_data_internal(*JJ, *II, *IF);
	      }
	    }
	    ++JJ; ++II;
	  }
	}
	++J;
      }
    }
  }

  void output_transient_state(Ioss::Region &output_region, RegionVector &part_mesh, double time,
			      const Map &local_node_map, SystemInterface &interface)
  {
    // Determine which state on each input mesh corresponds to 'time'
    std::vector<int> steps(part_mesh.size());
    for (size_t p=0; p < part_mesh.size(); p++) {
      size_t nts = part_mesh[p]->get_property("state_count").get_int();
      steps[p] = 0;
      for (size_t i=0; i < nts; i++) {
	if (approx_equal(part_mesh[p]->get_state_time(i+1), time)) {
	  steps[p] = i+1;
	  part_mesh[p]->begin_state(steps[p]);
	  break;
	}
      }
    }

    output_globals(output_region, part_mesh, time, steps);
    output_nodal(output_region, part_mesh, time, steps, local_node_map, interface);
    output_element(output_region, part_mesh, time, steps);
    output_nodal_nodeset_fields(output_region, part_mesh, time, steps, interface);
    if (!interface.omit_nodesets()) {
      output_nset(output_region, part_mesh, time, steps);
    }
    if (!interface.omit_sidesets()) {
      output_sset(output_region, part_mesh, time, steps);
    }

    for (size_t p=0; p < part_mesh.size(); p++) {
      if (steps[p] != 0)
	part_mesh[p]->end_state(steps[p]);
    }
  }

  void transfer_field_data(Ioss::GroupingEntity *ige,
                           Ioss::GroupingEntity *oge,
                           Ioss::Field::RoleType role,
                           const std::string &prefix,
                           bool transfer_connectivity)
  {
    // Iterate through the TRANSIENT-role fields of the input
    // database and transfer to output database.
    Ioss::NameList state_fields;
    Ioss::NameList::const_iterator IF;
    ige->field_describe(role, &state_fields);

    // Complication here is that if the 'role' is 'Ioss::Field::MESH',
    // then the 'ids' field must be transferred first...
    if (role == Ioss::Field::MESH) {
      for (IF = state_fields.begin(); IF != state_fields.end(); ++IF) {
        std::string field_name = *IF;
        assert(oge->field_exists(field_name));
        if (field_name == "ids") {
          transfer_field_data_internal(ige, oge, field_name);
          break;
        }
      }
    }

    for (IF = state_fields.begin(); IF != state_fields.end(); ++IF) {
      std::string field_name = *IF;
      // All of the 'Ioss::EntityBlock' derived classes have a
      // 'connectivity' field, but it is only interesting on the
      // Ioss::ElementBlock class. On the other classes, it just
      // generates overhead...
      if (!transfer_connectivity && field_name == "connectivity")
        continue;


      if (field_name != "ids" &&
          (prefix.length() == 0 || std::strncmp(prefix.c_str(), field_name.c_str(), prefix.length()) == 0)) {
        assert(oge->field_exists(field_name));
        transfer_field_data_internal(ige, oge, field_name);
      }
    }
  }

  void transfer_field_data_internal(Ioss::GroupingEntity *ige,
				    Ioss::GroupingEntity *oge,
				    const std::string &field_name)
  {

    size_t isize = ige->get_field(field_name).get_size();
    assert (isize == oge->get_field(field_name).get_size());

    std::vector<double> data;
    ige->get_field_data(field_name, data);
    oge->put_field_data(field_name, data);
  }

  void define_global_fields(Ioss::Region &output_region, RegionVector &part_mesh,
			    const StringIdVector &variable_list)
  {
    if (!variable_list.empty() && variable_list[0].first == "none")
      return;
    size_t part_count = part_mesh.size();
    for (size_t p=0; p<part_count; p++) {
      Ioss::NameList fields;
      part_mesh[p]->field_describe(Ioss::Field::TRANSIENT, &fields);
      Ioss::NameList::const_iterator IF;
      for (IF = fields.begin(); IF != fields.end(); ++IF) {
	if (valid_variable(*IF, 0, variable_list)) {
	  Ioss::Field field = part_mesh[p]->get_field(*IF);
	  output_region.field_add(field);
	}
      }
    }
  }

  void define_nodal_fields(Ioss::Region &output_region, RegionVector &part_mesh,
			   const StringIdVector &variable_list, SystemInterface &interface)
  {
    if (!variable_list.empty() && variable_list[0].first == "none")
      return;
    Ioss::NodeBlock *onb = output_region.get_node_blocks()[0];
    SMART_ASSERT(onb != NULL);
    int node_count = onb->get_property("entity_count").get_int();
    size_t part_count = part_mesh.size();
    for (size_t p=0; p<part_count; p++) {
      if (!interface.convert_nodes_to_nodesets(p+1)) {
	Ioss::NodeBlock *nb = part_mesh[p]->get_node_blocks()[0];
	Ioss::NameList fields;
	SMART_ASSERT(nb != NULL);
	nb->field_describe(Ioss::Field::TRANSIENT, &fields);
	Ioss::NameList::const_iterator IF;
	for (IF = fields.begin(); IF != fields.end(); ++IF) {
	  if (valid_variable(*IF, 0, variable_list)) {
	    Ioss::Field field = nb->get_field(*IF);
	    field.reset_count(node_count);
	    onb->field_add(field);
	  }
	}
      }
    }
  }

  void define_element_fields(Ioss::Region &output_region, RegionVector &part_mesh,
			     const StringIdVector &variable_list)
  {
    // Element Block Fields...
    if (!variable_list.empty() && variable_list[0].first == "none")
      return;
    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      Ioss::ElementBlockContainer iebs = part_mesh[p]->get_element_blocks();
      Ioss::ElementBlockContainer::const_iterator J = iebs.begin();
      while (J != iebs.end()) {
	std::string name = part_mesh[p]->name() + "_" + (*J)->name();
	Ioss::ElementBlock *oeb = output_region.get_element_block(name);
	if (oeb != NULL) {
	  int id = oeb->get_property("id").get_int();
	  Ioss::NameList fields;
	  (*J)->field_describe(Ioss::Field::TRANSIENT, &fields);
	  Ioss::NameList::const_iterator IF;
	  for (IF = fields.begin(); IF != fields.end(); ++IF) {
	    if (valid_variable(*IF, id, variable_list)) {
	      Ioss::Field field = (*J)->get_field(*IF);
	      oeb->field_add(field);
	    }
	  }
	}
	++J;
      }
    }
  }

  void define_nset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
			  const StringIdVector &variable_list)
  {
    // Nodeset fields...
    if (!variable_list.empty() && variable_list[0].first == "none")
      return;

    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      Ioss::NodeSetContainer ins = part_mesh[p]->get_nodesets();
      Ioss::NodeSetContainer::const_iterator J = ins.begin();
      while (J != ins.end()) {
	if (!entity_is_omitted(*J)) {
	  std::string name = part_mesh[p]->name() + "_" + (*J)->name();
	  Ioss::NodeSet *ons = output_region.get_nodeset(name);
	  SMART_ASSERT(ons != NULL)(name);

	  int id = (*J)->get_property("id").get_int();
	  Ioss::NameList fields;
	  (*J)->field_describe(Ioss::Field::TRANSIENT, &fields);
	  Ioss::NameList::const_iterator IF;
	  for (IF = fields.begin(); IF != fields.end(); ++IF) {
	    if (valid_variable(*IF, id, variable_list)) {
	      Ioss::Field field = (*J)->get_field(*IF);
	      ons->field_add(field);
	    }
	  }
	}
	++J;
      }
    }
  }

  void define_sset_fields(Ioss::Region &output_region, RegionVector &part_mesh,
			  const StringIdVector &variable_list)
  {
    if (!variable_list.empty() && variable_list[0].first == "none")
      return;
    Ioss::SideSetContainer os = output_region.get_sidesets();
    Ioss::SideSetContainer::const_iterator I = os.begin();

    Ioss::SideBlockContainer out_eb;
    // Put all output side blocks in the same list...
    while (I != os.end()) {
      Ioss::SideBlockContainer obs = (*I++)->get_side_blocks();
      std::copy(obs.begin(), obs.end(), std::back_inserter(out_eb));
    }
    
    // Assuming (with checks) that the output side blocks will be
    // iterated in same order as input side blocks...
    Ioss::SideBlockContainer::const_iterator II = out_eb.begin();

    size_t part_count = part_mesh.size();
    for (size_t p = 0; p < part_count; p++) {
      Ioss::SideSetContainer is = part_mesh[p]->get_sidesets();
      Ioss::SideSetContainer::const_iterator J = is.begin();
      while (J != is.end()) {
	if (!entity_is_omitted(*J)) {
	  int id = (*J)->get_property("id").get_int();
	  Ioss::SideBlockContainer ebs = (*J)->get_side_blocks();
	  Ioss::SideBlockContainer::const_iterator JJ = ebs.begin();

	  while (JJ != ebs.end()) {
	    SMART_ASSERT(part_mesh[p]->name() + "_" + (*JJ)->name() == (*II)->name());
	    Ioss::NameList fields;
	    (*JJ)->field_describe(Ioss::Field::TRANSIENT, &fields);
	    Ioss::NameList::const_iterator IF;
	    for (IF = fields.begin(); IF != fields.end(); ++IF) {
	      if (valid_variable(*IF, id, variable_list)) {
		Ioss::Field field = (*JJ)->get_field(*IF);
		(*II)->field_add(field);
	      }
	    }
	    ++JJ; ++II;
	  }
	}
	++J;
      }
    }
  }

  void transfer_fields(Ioss::GroupingEntity *ige,
                       Ioss::GroupingEntity *oge,
                       Ioss::Field::RoleType role,
                       const std::string &prefix)
  {
    // Check for transient fields...
    Ioss::NameList fields;
    ige->field_describe(role, &fields);

    // Iterate through results fields and transfer to output
    // database...  If a prefix is specified, only transfer fields
    // whose names begin with the prefix
    Ioss::NameList::const_iterator IF;
    for (IF = fields.begin(); IF != fields.end(); ++IF) {
      std::string field_name = *IF;
      if (field_name != "ids" && !oge->field_exists(field_name) &&
          (prefix.length() == 0 || std::strncmp(prefix.c_str(), field_name.c_str(), prefix.length()) == 0)) {
        // If the field does not already exist, add it to the output node block
        Ioss::Field field = ige->get_field(field_name);
        oge->field_add(field);
      }
    }
  }

  bool valid_variable(const std::string variable, int id, const StringIdVector &variable_list)
  {
    if (variable_list.empty() || variable_list[0].first == "all")
      return true;
    if (variable_list[0].first == "none")
      return false;

    StringIdVector::const_iterator IF;
    for (IF = variable_list.begin(); IF != variable_list.end(); ++IF) {
      if ((*IF).first == variable) {
	if (id == 0 || id == (*IF).second || (*IF).second == 0) {
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
	  Ioss::NodeSetContainer nodesets = part_mesh[p]->get_nodesets();
	  Ioss::NodeSetContainer::const_iterator I;
	  for (I=nodesets.begin(); I != nodesets.end(); ++I) {
	    (*I)->property_add(Ioss::Property(std::string("omitted"), 1));
	  }
	} else {
	  for (size_t nset = 0; nset < omit[p].size(); nset++) {
	    Ioss::NodeSet *ns = part_mesh[p]->get_nodeset(omit[p][nset]);
	    if (ns != NULL) {
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
	  Ioss::SideSetContainer sidesets = part_mesh[p]->get_sidesets();
	  Ioss::SideSetContainer::const_iterator I;
	  for (I=sidesets.begin(); I != sidesets.end(); ++I) {
	    (*I)->property_add(Ioss::Property(std::string("omitted"), 1));
	  }
	} else {
	  for (size_t sset = 0; sset < omit[p].size(); sset++) {
	    Ioss::SideSet *ss = part_mesh[p]->get_sideset(omit[p][sset]);
	    if (ss != NULL) {
	      ss->property_add(Ioss::Property(std::string("omitted"), 1));
	    }
	  }
	}
      }
    }    
  }
}
