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

#include <im_exodusII.h>
#include <im_ne_nemesisI.h>
#include <create_inline_mesh.h>

#include <pamgen/Iopg_DatabaseIO.h>

#include <Ioss_CodeTypes.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <iterator>
#include <time.h>
#include <float.h>
#include <limits.h>

namespace Iopg {
  typedef std::set<std::string> SideSetSet;
  typedef std::map<std::string, const std::string, std::less<const std::string> > SideSetMap;

  struct TopologyMapCompare {
    bool operator() (const std::pair<std::string, const Ioss::ElementTopology*> &lhs,
                     const std::pair<std::string, const Ioss::ElementTopology*> &rhs) const
    {
      assert(lhs.second != NULL);
      assert(rhs.second != NULL);
      return lhs.first < rhs.first || (!(rhs.first < lhs.first) && lhs.second->name() < rhs.second->name());
    }
  };

  typedef std::map<std::pair<std::string, const Ioss::ElementTopology*>, int, TopologyMapCompare > TopologyMap;
  typedef TopologyMap::value_type TopoMapPair;
}

namespace {
  // Output a message that the operation is unsupported and die...
  void unsupported(const char *operation)
  {
    std::cerr << "ERROR: Unsupported functionality called: " << operation << std::endl;
    std::abort();
  }
  
  const size_t max_string_length = MAX_STR_LENGTH;
  const size_t max_line_length   = MAX_LINE_LENGTH;
  const std::string SCALAR()     {return std::string("scalar");}
  const std::string VECTOR3D()   {return std::string("vector_3d");}
  const std::string SYM_TENSOR() {return std::string("sym_tensor_33");}

  void separate_surface_element_sides(Ioss::IntVector &element,
				      Ioss::IntVector &sides,
				      Ioss::Region *region,
				      Iopg::TopologyMap &topo_map,
				      Iopg::TopologyMap &side_map,
				      Ioss::SurfaceSplitType split_type);
    
  int get_file_pointer()
  {
    return 0;
  }

  const char *Version() {return "Iopg_DatabaseIO.C 2010/09/22";}

  void pamgen_error(int exoid, int lineno, int /* processor */) {
    std::ostringstream errmsg;

    errmsg << "Pamgen error at line " << lineno
	   << " in file '" << Version()
	   << "' Please report to gdsjaar@sandia.gov if you need help.";

    IOSS_ERROR(errmsg);
  }
}

namespace Iopg {
  // ========================================================================
  const IOFactory* IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory()
    : Ioss::IOFactory("pamgen")
  {}

  Ioss::DatabaseIO* IOFactory::make_IO(const std::string& filename,
				       Ioss::DatabaseUsage db_usage,
				       MPI_Comm communicator) const
  { return new DatabaseIO(NULL, filename, db_usage, communicator); }

  // ========================================================================
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
			 Ioss::DatabaseUsage db_usage,
			 MPI_Comm communicator) :
    Ioss::DatabaseIO(region, filename, db_usage, communicator), 
    spatialDimension(3), nodeCount(0),
    elementCount(0), nodeBlockCount(0),
    elementBlockCount(0), nodesetCount(0), sidesetCount(0),
    commsetNodeCount(0), commsetElemCount(0),
    sequentialNG2L(false), sequentialEG2L(false)
  {
    if (is_input()) {
      dbState = Ioss::STATE_UNKNOWN;
    } else {
      std::ostringstream errmsg;
      errmsg << "Pamgen mesh option is only valid for input mesh.";
      IOSS_ERROR(errmsg);
    }
  }

  DatabaseIO::~DatabaseIO()
  {
    Delete_Pamgen_Mesh();
  }

  void DatabaseIO::read_meta_data()
  {
    // The file for pamgen contains the mesh description.
    // The Iopg routine is expecting the mesh description to be a
    // single string...

    // Read the data and convert it to a string which is then passed
    // in to the Create_Pamgen_Mesh routine.
    std::ifstream input(get_filename().c_str());
    if (!input) {
      std::ostringstream errmsg;
      errmsg << "Error opening file '" << get_filename() << "'.";
      IOSS_ERROR(errmsg);
    }

    std::string mesh_description;
    std::string line;
    int dimension = 3;
    while (std::getline(input, line)) {
      if (line.empty() || line[0] == '$' || line[0] == '#') {
	continue;
      } else {
	bool twod = line.find('2') != std::string::npos;
	bool trid = line.find('3') != std::string::npos;
	if (twod && trid) {
	  std::ostringstream errmsg;
	  errmsg << "PAMGEN: Cannot determine spatial dimension; found both 2 and 3 on first line of file.";
	  IOSS_ERROR(errmsg);
	} else if (twod) {
	  dimension = 2;
	} else if (trid) {
	  dimension = 3;
	}
	break;
      }
    }

    mesh_description += "mesh\n";
    while (std::getline(input, line)) {
      if (line.empty() || line[0] == '$' || line[0] == '#') {
	;
      } else {
	mesh_description += line;
	mesh_description += "\n";
      }
    }

    int retval = ERROR_FREE_CREATION;
    bool error_detected = false;
    try {
      retval = Create_Pamgen_Mesh(mesh_description.c_str(), dimension,
				  util().parallel_rank(), util().parallel_size(), INT_MAX);
    }
    catch (const std::exception &x) {
      error_detected = true;
    }

    if (error_detected || retval != ERROR_FREE_CREATION) {
      char *error = NULL;
      if (retval == ERROR_CREATING_MS) {
	int err_size = getPamgenErrorStreamSize();
	error = new char[err_size+1];
	getPamgenErrorStream(error);
      } else {
	int err_size = getPamgenEchoStreamSize();
	error = new char[err_size+1];
	getPamgenEchoStream(error);
      }
      // handle error
      std::ostringstream errmsg;
      errmsg << "Pamgen mesh generation failed:\n";
      errmsg << error;
      delete [] error;
      IOSS_ERROR(errmsg);
    }

    read_region();
    read_communication_metadata();

    get_elemblocks();
    check_side_topology();

    get_nodeblocks();
    get_sidesets();

    get_nodesets();
    get_commsets();
  }

  void DatabaseIO::read_region()
  {
    // Add properties and fields to the 'owning' region.
    // Also defines member variables of this class...

    char dbtitle[max_line_length+1];
    std::memset(dbtitle, 0, max_line_length+1);

    int error = im_ex_get_init(get_file_pointer(), dbtitle, &spatialDimension,
			       &nodeCount, &elementCount, &elementBlockCount,
			       &nodesetCount, &sidesetCount);
    if (error < 0)
      pamgen_error(get_file_pointer(), __LINE__, myProcessor);

    nodeBlockCount = 1;

    if (nodeCount == 0) {
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      IOSS_WARNING << "No nodes were found in the model, file '" << decoded_filename << "'";
    } else if (nodeCount < 0) {
      // NOTE: Code will not continue past this call...
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      std::ostringstream errmsg;
      errmsg << "Negative node count was found in the model\n"
	     << "       File: '" << decoded_filename << "'.\n";
      IOSS_ERROR(errmsg);
    }

    if (elementCount == 0) {
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      IOSS_WARNING << "No elements were found in the model, file: '" << decoded_filename << "'";
    }

    if (elementCount < 0) {
      // NOTE: Code will not continue past this call...
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      std::ostringstream errmsg;
      errmsg << "Negative element count was found in the model, file: '"
	     << decoded_filename << "'";
      IOSS_ERROR(errmsg);
    }

    if (elementCount > 0 && elementBlockCount <= 0) {
      // NOTE: Code will not continue past this call...
      std::string decoded_filename = util().decode_filename(get_filename(), isParallel);
      std::ostringstream errmsg;
      errmsg << "No element blocks were found in the model, file: '" << decoded_filename << "'";
      IOSS_ERROR(errmsg);
    }

    Ioss::Region *this_region = get_region();
    this_region->property_add(Ioss::Property(std::string("title"), dbtitle));
    this_region->property_add(Ioss::Property(std::string("spatial_dimension"),
					     spatialDimension));
  }

  void DatabaseIO::read_communication_metadata()
  {
    // Sierra does not use/need the element communication map data,
    // but the Nemesis api call to get the nodal communication map
    // data also gets the element communication data. So,...,we need to
    // do some redundant reads to get back the element comm data
    // and also some extra memory allocations to give them a place
    // to go.

    // Check that file is nemesis.
    int num_proc;          // Number of processors file was decomposed for
    int num_proc_in_file;  // Number of processors this file has info for
    char file_type[2];         // "s" for scalar, "p" for parallel
    // Sierra is set up to handle 1 processor per file; parallel file

    int error = im_ne_get_init_info(get_file_pointer(),
				    &num_proc, &num_proc_in_file, &file_type[0]);
    if (error < 0) {
      // Not a nemesis file
      if (util().parallel_size() > 1) {
	std::ostringstream errmsg;
	errmsg << "Exodus file does not contain nemesis information.\n";
	IOSS_ERROR(errmsg);
      }
    } else if (num_proc != util().parallel_size() && util().parallel_size() > 1) {
      std::ostringstream errmsg;
      errmsg <<  "Exodus file was decomposed for " << num_proc
	     << " processors; application is currently being run on "
	     << util().parallel_size() << " processors";
      IOSS_ERROR(errmsg);

    } else if (num_proc_in_file != 1) {
      std::ostringstream errmsg;
      errmsg <<"Exodus file contains data for " << num_proc_in_file
	     << " processors; application requires 1 processor per file.";
      IOSS_ERROR(errmsg);

    } else if (file_type[0] != 'p') {
      std::ostringstream errmsg;
      errmsg << "Exodus file contains scalar nemesis data; application requires parallel nemesis data.";
      IOSS_ERROR(errmsg);
    }

    // Get global data (over all processors)
    int global_nodes       = nodeCount;
    int global_elements    = elementCount;
    int global_eblocks     = 0; // unused
    int global_nsets       = 0; // unused
    int global_ssets       = 0; // unused

    int num_external_nodes; // unused
    int num_elem_cmaps     = 0;
    int num_node_cmaps     = 0;
    int num_internal_nodes = nodeCount;
    int num_border_nodes   = 0;
    int num_internal_elems = elementCount;
    int num_border_elems   = 0;

    if (isParallel) {

      error = im_ne_get_loadbal_param(get_file_pointer(),
				      &num_internal_nodes,
				      &num_border_nodes,
				      &num_external_nodes,
				      &num_internal_elems,
				      &num_border_elems,
				      &num_node_cmaps,
				      &num_elem_cmaps,
				      myProcessor);
      if (error < 0)
	pamgen_error(get_file_pointer(), __LINE__, myProcessor);

      // A nemesis file typically separates nodes into multiple
      // communication sets by processor.  (each set specifies
      // nodes/elements that communicate with only a single processor).
      // For Sierra, we want a single node commun. map and a single
      // element commun. map specifying all communications so we combine
      // all sets into a single set.
      error = im_ne_get_init_global(get_file_pointer(), &global_nodes, &global_elements,
				    &global_eblocks, &global_nsets, &global_ssets);
      if (error < 0)
	pamgen_error(get_file_pointer(), __LINE__, myProcessor);
    }

    Ioss::Region *region = get_region();
    region->property_add(Ioss::Property("internal_node_count",
					num_internal_nodes));
    region->property_add(Ioss::Property("border_node_count",
					num_border_nodes));
    region->property_add(Ioss::Property("internal_element_count",
					num_internal_elems));
    region->property_add(Ioss::Property("border_element_count",
					num_border_elems));
    region->property_add(Ioss::Property("global_node_count",    global_nodes));
    region->property_add(Ioss::Property("global_element_count", global_elements));

    // Possibly, the following 4 fields should be nodesets and element
    // sets instead of fields on the region...
    region->field_add(Ioss::Field("internal_nodes", Ioss::Field::INTEGER,
				  SCALAR(),
				  Ioss::Field::COMMUNICATION,
				  num_internal_nodes));
    region->field_add(Ioss::Field("border_nodes", Ioss::Field::INTEGER,
				  SCALAR(),
				  Ioss::Field::COMMUNICATION,
				  num_border_nodes));
    region->field_add(Ioss::Field("internal_elements", Ioss::Field::INTEGER,
				  SCALAR(),
				  Ioss::Field::COMMUNICATION,
				  num_internal_elems));
    region->field_add(Ioss::Field("border_elements", Ioss::Field::INTEGER,
				  SCALAR(),
				  Ioss::Field::COMMUNICATION,
				  num_border_elems));

    assert(nodeCount    == num_internal_nodes + num_border_nodes);
    assert(elementCount == num_internal_elems + num_border_elems);
  }

  void DatabaseIO::get_nodeblocks() 
  {
    // For exodusII, there is only a single node block which contains
    // all of the nodes.
    // The default id assigned is '1' and the name is 'nodeblock_1'
    std::string block_name = "nodeblock_1";
    Ioss::NodeBlock *block = new Ioss::NodeBlock(this, block_name,
						 nodeCount, spatialDimension);
    block->property_add(Ioss::Property("id", 1));

    get_region()->add(block);
  }

  void DatabaseIO::get_elemblocks() 
  {
    // Attributes of an element block are:
    // -- id
    // -- name
    // -- element type
    // -- number of elements
    // -- number of attributes per element
    // -- number of nodes per element (derivable from type)
    // -- number of faces per element (derivable from type)
    // -- number of edges per element (derivable from type)

    // In a parallel execution, it is possible that an element block will have
    // no elements on a particular processor...

    // NOTE: This routine may be called multiple times on a single database.
    //       make sure it is not dependent on being called one time only...

    // Get exodusII element block metadata
    if (elementBlockCount == 0)
      return;

    Ioss::IntVector element_block_ids(elementBlockCount);

    int error = im_ex_get_elem_blk_ids(get_file_pointer(), &element_block_ids[0]);
    if (error < 0) {
      pamgen_error(get_file_pointer(), __LINE__, myProcessor);
    }

    size_t all_element_type_length = elementBlockCount * (max_string_length+1);
    std::vector<char> all_element_type(all_element_type_length);

    Ioss::IntVector attributes(elementBlockCount);
    Ioss::IntVector my_node_count(elementBlockCount);
    Ioss::IntVector local_element_count(elementBlockCount);
    Ioss::IntVector global_element_count(elementBlockCount);
    int iblk ;

    for ( iblk = 0; iblk < elementBlockCount; iblk++) {
      int id = element_block_ids[iblk];
      int number_elements;
      int nodes_per_element;
      int attributes_per_element;
      
      char * const element_type = &all_element_type[0] + iblk * (max_string_length+1);
      
      error = im_ex_get_elem_block(get_file_pointer(), id,
				   element_type,
				   &number_elements, &nodes_per_element,
				   &attributes_per_element);
      if (error < 0) {
	pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }
      
      local_element_count[iblk] = number_elements;
      
      if (number_elements == 0) {
	attributes[iblk] = 0;
	my_node_count[iblk] = 0;
      } else {
	attributes[iblk] = attributes_per_element;
	my_node_count[iblk] = nodes_per_element;
      }
    }

    // Determine global element count for each element block....
    // Can also get this from an nemesis call, but the data may not always be there in all cases.
    util().global_count(local_element_count, global_element_count);
    
    // The 'offset' is used to map an element location within an element
    // block to the element 'file descriptor'.  For example, the file
    // descriptor of the 37th element in the 4th block is calculated by:
    // file_descriptor = offset of block 4 + 37 This can also be used to
    // determine which element block an element with a file_descriptor
    // maps into. An particular element block contains all elements in
    // the range:
    //     offset < file_descriptor <= offset+number_elements_per_block
    int offset = 0;
    int used_blocks = 0;
    
    for ( iblk = 0; iblk < elementBlockCount; iblk++) {
      int id = element_block_ids[iblk];
      std::string alias = Ioss::Utils::encode_entity_name("block", id);
      char * const element_type = &all_element_type[0] + iblk * (max_string_length+1);

      Ioss::ElementBlock *block = NULL;
      std::string block_name = Ioss::Utils::encode_entity_name("block", id);

      Ioss::Utils::fixup_name(element_type); // Convert to lowercase; replace spaces with '_'
      std::string type = std::string(element_type);
      std::string save_type = type;

      // Fixup an exodusII kluge/ambiguity.
      // The element block type does not fully define the element. For
      // example, a block of type 'triangle' may have either 3 or 6
      // nodes.  To fix this, check the block type name and see if it
      // ends with a number.  If it does, assume it is OK; if not, append
      // the 'nodes_per_element'.
      if (!isdigit(*(type.rbegin()))) {
	if (my_node_count[iblk] > 1) {
	  type += Ioss::Utils::to_string(my_node_count[iblk]);
	}
      }

      block = new Ioss::ElementBlock(this, block_name, type,
				     local_element_count[iblk],
				     attributes[iblk]);

      block->property_add(Ioss::Property("id", id));
      
      // Maintain block order on output database...
      block->property_add(Ioss::Property("original_block_order", used_blocks++));
      
      if (block->get_property("topology_type").get_string() != save_type
	  && save_type != "null" && save_type != "") {
	// Maintain original element type on output database if possible.
	block->property_add(Ioss::Property("original_element_type", save_type));
      }
      
      block->property_add(Ioss::Property("global_entity_count", global_element_count[iblk]));
      
      offset += local_element_count[iblk];
      
      get_region()->add(block);
      if (block_name != alias) {
	get_region()->add_alias(block_name, alias);
      }
    }
    assert(elementCount == offset);
  }

  void DatabaseIO::get_nodesets()
  {
    // Attributes of a nodeset are:
    // -- id
    // -- name
    // -- number of nodes
    // -- number of distribution factors (see next comment)
    // ----the #distribution factors should equal #nodes or 0, any
    //     other value does not make sense. If it is 0, then a substitute
    //     list will be created returning 1.0 for the factor

    // In a parallel execution, it is possible that a nodeset will have
    // no nodes or distribution factors on a particular processor...

    // Get exodusII nodeset metadata
    if (nodesetCount > 0) {
      Ioss::IntVector nodeset_ids(nodesetCount);
      int error = im_ex_get_node_set_ids(get_file_pointer(), &nodeset_ids[0]);
      if (error < 0) {
	pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }

      for (int ins = 0; ins < nodesetCount; ins++) {
	int id = nodeset_ids[ins];
	int number_nodes;
	int number_distribution_factors;

	error = im_ex_get_node_set_param(get_file_pointer(), id,
					 &number_nodes, &number_distribution_factors);
	if (error < 0) {
	  pamgen_error(get_file_pointer(), __LINE__, myProcessor);
	}

	std::string nodeset_name = Ioss::Utils::encode_entity_name("nodelist", id);
	Ioss::NodeSet *nodeset = new Ioss::NodeSet(this, nodeset_name, number_nodes);
	nodeset->property_add(Ioss::Property("id", id));
	get_region()->add(nodeset);

	get_region()->add_alias(nodeset_name, Ioss::Utils::encode_entity_name("nodelist", id));
      }

    }
  }

  void DatabaseIO::get_commsets()
  {
    // Attributes of a commset are:
    // -- id (property)
    // -- name (property)
    // -- number of node--CPU pairs (field)

    // In a parallel execution, it is possible that a commset will have
    // no nodes on a particular processor...

    // If this is a serial execution, there will be no communication
    // nodesets, just return an empty container.

    if (isParallel) {
      // This is a parallel run. There should be communications data
      // Get nemesis commset metadata
      int my_node_count = 0;
      int elem_count = 0;

      // NOTE: It is possible for a parallel run to have no
      // communications maps if the decomposition occurs along contact
      // surfaces.  In this case, we create empty node and element
      // communication maps.
      if (commsetNodeCount > 0 || commsetElemCount > 0) {
	if (commsetNodeCount > 0) {
	  nodeCmapIds      = new int[commsetNodeCount];
	  nodeCmapNodeCnts = new int[commsetNodeCount];
	}
	if (commsetElemCount > 0) {
	  elemCmapIds      = new int[commsetElemCount];
	  elemCmapElemCnts = new int[commsetElemCount];
	}

	int error = im_ne_get_cmap_params(get_file_pointer(),
					  nodeCmapIds, nodeCmapNodeCnts,
					  elemCmapIds, elemCmapElemCnts,
					  myProcessor);
	if (error < 0)
	  pamgen_error(get_file_pointer(), __LINE__, myProcessor);

	// Count nodes, elements, and convert counts to offsets.
	//
	for (int ics = 0; ics < commsetNodeCount; ics++) {
	  my_node_count += nodeCmapNodeCnts[ics];
	}
	for (int ecs = 0; ecs < commsetElemCount; ecs++) {
	  elem_count += elemCmapElemCnts[ecs];
	}
      }
      // Create a single node commset and a single element commset
      Ioss::CommSet *commset = new Ioss::CommSet(this, "commset_node", "node",
						 my_node_count);
      commset->property_add(Ioss::Property("id", 1));
      get_region()->add(commset);

      commset = new Ioss::CommSet(this, "commset_side", "side",
				  elem_count);
      commset->property_add(Ioss::Property("id", 1));
      get_region()->add(commset);
    }
  }

  void DatabaseIO::get_sidesets()
  {
    // This function creates all sidesets (surfaces) for a
    // model.  Note that a side set contains 1 or more side
    // blocks which are homogenous (same topology). In serial execution,
    // this is fairly straightforward since there are no null sets and
    // we have all the information we need. (...except see below for
    // surface evolution).
    //
    // However, in a parallel execution, we have the possibility that a
    // side set will have no edges/faces or distribution factors on
    // a particular processor.  We then don't know the block topology of
    // the block(s) contained in this set. We could do some
    // communication and get a good idea of the topologies that are in
    // the set.
    // 
    // for now, If the database is three-dimensional, then sidesets are
    // assumed to translate into 'sidesets'; if two-dimensional, then
    // 'edgesets'

    if (sidesetCount > 0) {
      check_side_topology();

      // Get exodusII sideset metadata

      Ioss::IntVector side_set_ids(sidesetCount);
      int error = im_ex_get_side_set_ids(get_file_pointer(), &side_set_ids[0]);
      if (error < 0) {
	pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }

      for (int iss = 0; iss < sidesetCount; iss++) {
	int id = side_set_ids[iss];
	std::string sid = "";
	int number_sides;
	int number_distribution_factors;
	TopologyMap topo_map;
	TopologyMap side_map; // Used to determine side consistency

	Ioss::SurfaceSplitType split_type = splitType;
	std::string side_set_name;

	side_set_name = Ioss::Utils::encode_entity_name("surface", id);

	Ioss::SideSet *side_set = new Ioss::SideSet(this, side_set_name);
	get_region()->add(side_set);
	side_set->property_add(Ioss::Property("id", id));
	
	get_region()->add_alias(side_set_name, Ioss::Utils::encode_entity_name("surface", id));
	get_region()->add_alias(side_set_name, Ioss::Utils::encode_entity_name("sideset", id));

	//	  split_type = SPLIT_BY_ELEMENT_BLOCK;
	//	  split_type = SPLIT_BY_TOPOLOGIES;
	//	  split_type = SPLIT_BY_DONT_SPLIT;

	// Determine how many side blocks compose this side set.
	error = im_ex_get_side_set_param(get_file_pointer(), id,
					 &number_sides, &number_distribution_factors);
	if (error < 0) {
	  pamgen_error(get_file_pointer(), __LINE__, myProcessor);
	}

	Ioss::IntVector element(number_sides);
	Ioss::IntVector sides(number_sides);
	  
	int ierr = im_ex_get_side_set(get_file_pointer(), id, &element[0], &sides[0]);
	if (ierr < 0)
	  pamgen_error(get_file_pointer(), __LINE__, myProcessor);
	  
	if (split_type == Ioss::SPLIT_BY_TOPOLOGIES && sideTopology.size() == 1) {
	  // There is only one side type for all elements in the model
	  topo_map[std::make_pair(sideTopology[0].first->name(), sideTopology[0].second)] = number_sides;
	  
	} else if (split_type == Ioss::SPLIT_BY_DONT_SPLIT) {
	  const Ioss::ElementTopology *mixed_topo = Ioss::ElementTopology::factory("unknown");
	  topo_map[std::make_pair("unknown",mixed_topo)] = number_sides;
	  
	} else if (split_type == Ioss::SPLIT_BY_TOPOLOGIES) {
	  // There are multiple side types in the model.
	  // Iterate through the elements in the sideset, determine
	  // their parent element block using the blocks element
	  // topology and the side number, determine the side
	  // type.
	  
	  for (Ioss::TopoContainer::size_type i=0; i < sideTopology.size(); i++) {
	    topo_map[std::make_pair(sideTopology[i].first->name(), sideTopology[i].second)] = 0;
	    side_map[std::make_pair(sideTopology[i].first->name(), sideTopology[i].second)] = 0;
	  }
	  
	  separate_surface_element_sides(element, sides, get_region(), topo_map, side_map, split_type);
	}
	else if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
	  // There are multiple side types in the model.  Iterate
	  // through the elements in the sideset, determine their parent
	  // element block using blocks element topology and the side
	  // number, determine the side type.

	  // Seed the topo_map map with <block->name, face_topo>
	  // pairs so we are sure that all processors have the same
	  // starting topo_map (size and order).
	  Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
	  
	  for (int i=0; i < elementBlockCount; i++) {
	    Ioss::ElementBlock *block = element_blocks[i];
	    std::string name = block->name();
	    const Ioss::ElementTopology *common_ftopo = block->topology()->boundary_type(0);
	    if (common_ftopo != NULL) {
	      // All sides of this element block's topology have the same topology
	      topo_map[std::make_pair(name, common_ftopo)] = 0;
	      side_map[std::make_pair(name, common_ftopo)] = 0;
	    } else {
	      // The sides have different topology, iterate over
	      // them and create an entry for the unique face
	      // topology types
	      int par_dim = block->topology()->parametric_dimension();
	      if (par_dim == 2 || par_dim == 3) {
		int my_side_count = block->topology()->number_boundaries();
		for (int ii = 0; ii < my_side_count; ii++) {
		  const Ioss::ElementTopology *topo = block->topology()->boundary_type(ii+1);
		  topo_map[std::make_pair(name, topo)] = 0;
		  side_map[std::make_pair(name, topo)] = 0;
		}
	      }
	    }
	  }
	  separate_surface_element_sides(element, sides, get_region(), topo_map, side_map, split_type);
	}

	// End of first step in splitting.  Check among all processors
	// to see which potential splits have sides in them...
	Ioss::IntVector global_side_counts(topo_map.size());
	{
	  int i = 0;
	  {
	    TopologyMap::const_iterator I = topo_map.begin();
	    while (I != topo_map.end()) {
	      global_side_counts[i++] = (*I).second;
	      ++I;
	    }
	  }

	  // If splitting by element block, also sync the side_map
	  // information which specifies whether the sideset has
	  // consistent sides for all elements. Only really used for
	  // shells, but easier to just set the value on all surfaces
	  // in the element block split case.
	  if (side_map.size() == topo_map.size()) {
	    global_side_counts.resize(topo_map.size() + side_map.size());

	    {
	      TopologyMap::const_iterator I = side_map.begin();
	      while (I != side_map.end()) {
		global_side_counts[i++] = (*I).second;
		++I;
	      }
	    }
	  }

	  // See if any processor has non-zero count for the topo_map counts
	  // For the side_map, need the max value.
	  util().global_array_minmax(&global_side_counts[0], global_side_counts.size(),
				     Ioss::ParallelUtils::DO_MAX);
	}


	// Create Side Blocks

	int i = 0;
	TopologyMap::const_iterator I = topo_map.begin();
	assert(I != topo_map.end());
	while (I != topo_map.end()) {
	  if (global_side_counts[i++] > 0) {
	    const std::string topo_or_block_name   = (*I).first.first;
	    const Ioss::ElementTopology *side_topo = (*I).first.second;
	    assert(side_topo != NULL);
#if 0
	    if (side_topo->parametric_dimension() == topology_dimension-1 ||
		split_type == Ioss::SPLIT_BY_DONT_SPLIT ) {
#else
	    if (true) {
#endif
	      int side_count = (*I).second;

	      std::string side_block_name = "surface_" + topo_or_block_name + "_" + side_topo->name();
	      if (side_set_name == "universal_sideset") {
		side_block_name = side_set_name;
	      } else {
		if (sid == "")
		  side_block_name = Ioss::Utils::encode_entity_name(side_block_name, id);
		else {
		  side_block_name += "_";
		  side_block_name += sid;
		}
	      }

	      Ioss::ElementBlock* block = NULL;
	      // Need to get elem_topo....
	      const Ioss::ElementTopology *elem_topo = NULL;
	      if (split_type == Ioss::SPLIT_BY_TOPOLOGIES) {
		elem_topo = Ioss::ElementTopology::factory(topo_or_block_name);
	      }
	      else if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
		block = get_region()->get_element_block(topo_or_block_name);
		if (block == NULL) {
		  std::ostringstream errmsg;
		  std::cerr << "INTERNAL ERROR: Could not find element block '" << topo_or_block_name 
			    << "' Something is wrong in the Ioex::DatabaseIO class. Please report.\n";
		  IOSS_ERROR(errmsg);
		}
		elem_topo = block->topology();
	      }
	      if (split_type == Ioss::SPLIT_BY_DONT_SPLIT) {
		// Most likely this is "unknown", but can be a true
		// topology if there is only a single element block in
		// the model.
		elem_topo = Ioss::ElementTopology::factory(topo_or_block_name);
	      }
	      assert(elem_topo != NULL);

	      Ioss::SideBlock *side_block = new Ioss::SideBlock(this, side_block_name,
								side_topo->name(),
								elem_topo->name(),
								side_count);
	      side_set->add(side_block);

	      // Note that all sideblocks within a specific
	      // sideset might have the same id.
	      assert(side_block != NULL);
	      side_block->property_add(Ioss::Property("id", id));

	      // If splitting by element block, need to set the
	      // element block member on this side block.
	      if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
		side_block->set_parent_element_block(block);

	      }

	      // If we calculated whether the element side is
	      // consistent for all faces in this block, then
	      // tell the block which side it is, or that they are
	      // inconsistent. If it wasn't calculated above, then it
	      // will be calculated on the fly when/if requested.
	      // This is to avoid reading the sideset bulk data in
	      // cases where we don't need to read it, but if we are
	      // already reading it (to split the sidesets), then use
	      // the data when we have it.
	      if (side_map.size() > 0) {
		// Set a property indicating which element side
		// (1-based) all faces in this block are applied to.
		// If they are not all assigned to the same element
		// side, indicate this with a side equal to 0.
		//
		// (note: 'i' has already been incremented earlier in
		// the loop.  We need previous value here...)
		int side = global_side_counts[i-1+topo_map.size()];
		if (side == 999) side = 0;
		assert(side <= elem_topo->number_boundaries());
		side_block->set_consistent_side_number(side);
	      }

	      // Add an alias...
	      get_region()->add_alias(side_block);

	      if (split_type != Ioss::SPLIT_BY_DONT_SPLIT
		  && (number_distribution_factors > 0 || isParallel)
		  && side_set_name != "universal_sideset") {
		std::string storage = "Real[";
		storage += Ioss::Utils::to_string(side_topo->number_nodes());
		storage += "]";
		side_block->field_add(Ioss::Field("distribution_factors",
						Ioss::Field::REAL, storage,
						Ioss::Field::MESH, side_count));
	      }

	      if (side_set_name == "universal_sideset") {
		side_block->field_add(Ioss::Field("side_ids",
						Ioss::Field::INTEGER, "scalar",
						Ioss::Field::MESH, side_count));
	      }
	    }
	  }
	  ++I;
	}
      }
    }
  }

  bool DatabaseIO::begin(Ioss::State /* state */)
  {
    return true;
  }

  bool   DatabaseIO::end(Ioss::State /* state */)
  {
    return true;
  }

  bool DatabaseIO::begin_state(Ioss::Region *region, int /* state */, double time )
  {
    return true;
  }

  bool   DatabaseIO::end_state(Ioss::Region */* region */, int /* state */, double /* time */)
  {
    return true;
  }

  int DatabaseIO::get_field_internal(const Ioss::NodeBlock* nb,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {

      Ioss::Field::RoleType role = field.get_role();
      if (role == Ioss::Field::MESH) {
	if (field.get_name() == "mesh_model_coordinates") {
	  // Data required by upper classes store x0, y0, z0, ... xn,
	  // yn, zn. Data stored in exodusII file is x0, ..., xn, y0,
	  // ..., yn, z0, ..., zn so we have to allocate some scratch
	  // memory to read in the data and then map into supplied
	  // 'data'
	  std::vector<double> x(num_to_get);
	  std::vector<double> y(num_to_get);
	  std::vector<double> z;
	  if (spatialDimension == 3)
	    z.resize(num_to_get);

	  // Cast 'data' to correct size -- double
	  double *rdata = static_cast<double*>(data);

	  int ierr = im_ex_get_coord(get_file_pointer(), &x[0], &y[0], &z[0]);
	  if (ierr < 0)
	    pamgen_error(get_file_pointer(), __LINE__, myProcessor);

	  int index = 0;
	  for (size_t i=0; i < num_to_get; i++) {
	    rdata[index++] = x[i];
	    rdata[index++] = y[i];
	    if (spatialDimension == 3)
	      rdata[index++] = z[i];
	  }
	}

	else if (field.get_name() == "ids") {
	  // Map the local ids in this node block
	  // (1...node_count) to global node ids.
	  const Ioss::MapContainer &map = get_node_map();
	  int *ids = static_cast<int*>(data);

	  for (size_t i=0; i < num_to_get; i++) {
	    ids[i] = map[1 + i];
	  }
	}
	else if (field.get_name() == "connectivity") {
	  // Do nothing, just handles an idiosyncracy of the GroupingEntity
	}
	else {
	  num_to_get = Ioss::Utils::field_warning(nb, field, "input");
	}
      }
    }
    return num_to_get;
  }

  int DatabaseIO::get_field_internal(const Ioss::Region* /* reg */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  int DatabaseIO::get_field_internal(const Ioss::ElementBlock* eb, const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {
      
      int ierr = 0;
      int id = eb->get_property("id").get_int();
      int my_element_count = eb->get_property("entity_count").get_int();
      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
	// Handle the MESH fields required for an ExodusII file model.
	// (The 'genesis' portion)

	if (field.get_name() == "connectivity") {
	  int element_nodes = eb->get_property("topology_node_count").get_int();
	  assert(field.raw_storage()->component_count() == element_nodes);

	  // The connectivity is stored in a 1D array.
	  // The element_node index varies fastet
	  if (my_element_count > 0) {
	    ierr = im_ex_get_elem_conn(get_file_pointer(), id, static_cast<int*>(data));
	    if (ierr < 0)
	      pamgen_error(get_file_pointer(), __LINE__, myProcessor);

	    // Now, map the nodes in the connectivity from local to global ids
	    const Ioss::MapContainer &map = get_node_map();
	    if (!Ioss::Map::is_sequential(map)) {
	      int *connect = static_cast<int*>(data);
	      for (size_t i=0; i < num_to_get*element_nodes; i++)
		connect[i] = map[connect[i]];
	    }
	  }
	}
	else if (field.get_name() == "ids") {
	  // Map the local ids in this element block
	  // (eb_offset+1...eb_offset+1+my_element_count) to global element ids.
	  const Ioss::MapContainer &map = get_element_map();
	  int eb_offset = eb->get_offset();
	  int *ids = static_cast<int*>(data);

	  for (size_t i=0; i < num_to_get; i++) {
	    ids[i] = map[eb_offset + 1 + i];
	  }
	}
	else {
	  num_to_get = Ioss::Utils::field_warning(eb, field, "input");
	}
      }

      else if (role == Ioss::Field::ATTRIBUTE) {
	// ERROR: Pamgen doesn't support elements with attributes...
      } else if (role == Ioss::Field::TRANSIENT) {
	// ERROR: Pamgen doesn't support transient data...
      } else if (role == Ioss::Field::REDUCTION) {
	num_to_get = Ioss::Utils::field_warning(eb, field, "input reduction");
      }
    }
    return num_to_get;
  }

  int DatabaseIO::get_field_internal(const Ioss::CommSet* cs,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    {
      Ioss::SerializeIO	serializeIO__(this);

      size_t num_to_get = field.verify(data_size);

      if (num_to_get > 0) {
	int entity_count = cs->get_property("entity_count").get_int();

	// Return the <entity (node or side), processor> pair
	if (field.get_name() == "entity_processor") {

	  // Check type -- node or side
	  std::string type = cs->get_property("entity_type").get_string();

	  // Allocate temporary storage space
	  Ioss::IntVector entities(num_to_get);
	  Ioss::IntVector procs(num_to_get);

	  if (type == "node") {
	    int cm_offset = 0;
	    int i;
	    for (i=0; i < commsetNodeCount; i++) {
	      int ierr = im_ne_get_node_cmap(get_file_pointer(), nodeCmapIds[i],
					     &entities[cm_offset], &procs[cm_offset],
					     myProcessor);
	      if (ierr < 0)
		pamgen_error(get_file_pointer(), __LINE__, myProcessor);
	      cm_offset += nodeCmapNodeCnts[i];
	    }
	    assert(cm_offset == entity_count);

	    // Convert local node id to global node id and store in 'data'
	    int* entity_proc = static_cast<int*>(data);
	    const Ioss::MapContainer &map = get_node_map();

	    int j=0;
	    for (i=0; i < entity_count; i++) {
	      int local_id = entities[i];
	      entity_proc[j++] = map[local_id];
	      entity_proc[j++] = procs[i];
	    }
	  } else if (type == "side") {
	    Ioss::IntVector sides(entity_count);
	    int cm_offset = 0;
	    int i;
	    for (i=0; i < commsetElemCount; i++) {
	      int ierr = im_ne_get_elem_cmap(get_file_pointer(), elemCmapIds[i],
					     &entities[cm_offset], &sides[cm_offset],
					     &procs[cm_offset], myProcessor);
	      if (ierr < 0)
		pamgen_error(get_file_pointer(), __LINE__, myProcessor);
	      cm_offset += elemCmapElemCnts[i];
	    }
	    assert(cm_offset == entity_count);

	    int* entity_proc = static_cast<int*>(data);
	    const Ioss::MapContainer &map = get_element_map();

	    int j=0;
	    for (i=0; i < entity_count; i++) {
	      int local_id = entities[i];
	      entity_proc[j++] = 10*map[local_id]+sides[i];
	      entity_proc[j++] = procs[i];
	    }
	  } else {
	    std::ostringstream errmsg;
	    errmsg << "Invalid commset type " << type;
	    IOSS_ERROR(errmsg);
	  }

	} else {
	  num_to_get = Ioss::Utils::field_warning(cs, field, "input");
	}
      }
      return num_to_get;
    }
  }

  int DatabaseIO::get_field_internal(const Ioss::SideBlock* fb,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {

      int id = fb->get_property("id").get_int();
      size_t entity_count = fb->get_property("entity_count").get_int();
      if (num_to_get != entity_count) {
	std::ostringstream errmsg;
	errmsg << "Partial field input not yet implemented for side blocks";
	IOSS_ERROR(errmsg);
      }

      int number_sides;
      int number_distribution_factors;
      int ierr = im_ex_get_side_set_param(get_file_pointer(), id,
					  &number_sides, &number_distribution_factors);
      if (ierr < 0)
	pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      
      Ioss::Field::RoleType role = field.get_role();
      if (role == Ioss::Field::MESH) {


	// In exodusII, we may have split the sideset into multiple
	// side blocks if there are multiple side  topologies in the
	// sideset.  Because of this, the passed in 'data' may not be
	// large enough to hold the data residing in the sideset and we
	// may need to allocate a temporary array...  This can be checked
	// by comparing the size of the sideset with the 'side_count' of
	// the side block.

	if (field.get_name() == "side_ids") {
	}

	else if (field.get_name() == "ids") {
	  // In exodusII, the 'side set' is stored as a sideset.  A
	  // sideset has a list of elements and a corresponding local
	  // element side (1-based) The side id is: side_id =
	  // 10*element_id + local_side_number This assumes that all
	  // sides in a sideset are boundary sides.  Since we
	  // only have a single array, we need to allocate an extra array
	  // to store all of the data.  Note also that the element_id is
	  // the global id but only the local id is stored so we need to
	  // map from local_to_global prior to generating the side  id...

	  // Get the element number map (1-based)...
	  const Ioss::MapContainer &map = get_element_map();

	  // Allocate space for local side number, use 'data' as temporary
	  // storage for element numbers and overwrite with the side
	  // numbers.
	  Ioss::IntVector sides;
	  int *element = NULL;
	  int *ids = static_cast<int*>(data);
	  if (number_sides == static_cast<int>(entity_count)) {
	    // Only 1 side block in this sideset
	    sides.resize(entity_count);
	    element = static_cast<int*>(data);
	  } else {
	    // Multiple side blocks in sideset. Need to allocate memory to
	    // store sideset
	    sides.resize(number_sides);
	    element = new int[number_sides];
	  }

	  ierr = im_ex_get_side_set(get_file_pointer(), id, element, &sides[0]);
	  if (ierr < 0)
	    pamgen_error(get_file_pointer(), __LINE__, myProcessor);

	  if (number_sides == static_cast<int>(entity_count)) {
	    for (size_t iel = 0; iel < entity_count; iel++) {
	      int new_id = 10*map[element[iel]] + sides[iel];
	      ids[iel] = new_id;
	    }
	  } else {
	    Ioss::IntVector is_valid_side(number_sides);
	    Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, element, &sides[0],
							number_sides, get_region());
	    size_t ieb = 0;
	    for (int iel = 0; iel < number_sides; iel++) {
	      if (is_valid_side[iel] == 1) {
		int new_id = 10*map[element[iel]] + sides[iel];
		ids[ieb++] = new_id;
	      }
	    }
	    assert(ieb == entity_count);
	  }

	  if (number_sides != static_cast<int>(entity_count))
	    delete [] element;

	} else if (field.get_name() == "element_side") {
	  // In exodusII, the 'side set' is stored as a sideset.  A sideset
	  // has a list of elements and a corresponding local element side
	  // (1-based)

	  // Since we only have a single array, we need to allocate an extra
	  // array to store all of the data.  Note also that the element_id
	  // is the global id but only the local id is stored so we need to
	  // map from local_to_global prior to generating the side  id...

	  // Get the element number map (1-based)...
	  const Ioss::MapContainer &map = get_element_map();

	  // Allocate space for local side number and element numbers
	  // numbers.
	  int *element_side = static_cast<int*>(data);

	  Ioss::IntVector element(number_sides);
	  Ioss::IntVector sides(number_sides);

	  ierr = im_ex_get_side_set(get_file_pointer(), id, &element[0], &sides[0]);
	  if (ierr < 0)
	    pamgen_error(get_file_pointer(), __LINE__, myProcessor);

	  if (number_sides == static_cast<int>(entity_count)) {
	    size_t index = 0;
	    for (size_t iel = 0; iel < entity_count; iel++) {
	      element_side[index++] = map[element[iel]];
	      element_side[index++] = sides[iel];
	    }
	    assert(index/2 == entity_count);
	  } else {
	    Ioss::IntVector is_valid_side(number_sides);
	    Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, &element[0], &sides[0],
							number_sides, get_region());

	    size_t index = 0;
	    for (int iel = 0; iel < number_sides; iel++) {
	      if (is_valid_side[iel] == 1) {
		// This side  belongs in the side block
		element_side[index++] = map[element[iel]];
		element_side[index++] = sides[iel];
	      }
	    }
	    assert(index/2 == entity_count);
	  }

	} else if (field.get_name() == "connectivity") {
	  // The side  connectivity needs to be generated 'on-the-fly' from
	  // the element number and local side  of that element. A sideset
	  // can span multiple element blocks, and contain multiple side
	  // types; the side block contains side of similar topology.
	  ierr = get_side_connectivity(fb, id, entity_count, 
				       static_cast<int*>(data), data_size/sizeof(int));
	  if (ierr < 0)
	    pamgen_error(get_file_pointer(), __LINE__, myProcessor);

	} else if (field.get_name() == "distribution_factors") {
	  ierr = get_side_distributions(fb, id, entity_count, 
					static_cast<double*>(data), data_size/sizeof(double));
	  if (ierr < 0)
	    pamgen_error(get_file_pointer(), __LINE__, myProcessor);
	} else {
	  num_to_get = Ioss::Utils::field_warning(fb, field, "input");
	}
      } else if (role == Ioss::Field::TRANSIENT) {
	unsupported("side set transient fields");
      }
    }
    return num_to_get;
  }

  int DatabaseIO::get_field_internal(const Ioss::NodeSet* ns,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {

      int id = ns->get_property("id").get_int();
      Ioss::Field::RoleType role = field.get_role();
      if (role == Ioss::Field::MESH) {

	if (field.get_name() == "ids") {
	  int ierr = im_ex_get_node_set(get_file_pointer(), id, static_cast<int*>(data));
	  if (ierr < 0)
	    pamgen_error(get_file_pointer(), __LINE__, myProcessor);

	  // Convert the local node ids to global ids
	  const Ioss::MapContainer &map = get_node_map();
	  if (!Ioss::Map::is_sequential(map)) {
	    int *ids = static_cast<int*>(data);
	    for (size_t i=0; i < num_to_get; i++) {
	      ids[i] = map[ids[i]];
	    }
	  }
	} else if (field.get_name() == "distribution_factors") {
	  // Not supported by pamgen.  Fill data with '1'
	  double *rdata = static_cast<double*>(data);
	  for (size_t i=0; i < num_to_get; i++)
	    rdata[i] = 1.0;
	} else {
	  num_to_get = Ioss::Utils::field_warning(ns, field, "input");
	}
      } else if (role == Ioss::Field::TRANSIENT) {
	unsupported("transient nodeset fields");
      }
    }
    return num_to_get;
  }
    
  int DatabaseIO::get_field_internal(const Ioss::SideSet* /* fs */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::Region* /* region */, const Ioss::Field& /* field */,
				     void * /* data */, size_t /* data_size */) const
  {
    unsupported("output region field");
    return -1;
  }

  int DatabaseIO::put_field_internal(const Ioss::ElementBlock* /* eb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    unsupported("output element block field");
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::SideBlock* /* fb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    unsupported("output sideblock field");
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::NodeBlock* /* nb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    unsupported("output nodeblock field");
    return -1;
  }

  int DatabaseIO::put_field_internal(const Ioss::NodeSet* /* ns */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    unsupported("output nodeset field");
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::SideSet* /* fs */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    unsupported("output sideset field");
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::CommSet* /* cs */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    unsupported("output commset field");
    return -1;
  }

  const Ioss::MapContainer& DatabaseIO::get_node_map() const
  {
    // Allocate space for node number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (nodeMap.empty()) {
      nodeMap.resize(nodeCount+1);

      if (is_input()) {
	int error = im_ex_get_node_num_map(get_file_pointer(), &nodeMap[1]);
	if (error < 0) {
	  // Clear out the vector...
	  Ioss::MapContainer().swap(nodeMap);
	  pamgen_error(get_file_pointer(), __LINE__, myProcessor);
	}
	// Check for sequential node map.
	// If not, build the reverse G2L node map...
	sequentialNG2L = true;
	for (int i=1; i < nodeCount+1; i++) {
	  if (i != nodeMap[i]) {
	    sequentialNG2L = false;
	    nodeMap[0] = 1;
	    break;
	  }
	}

	if (!sequentialNG2L) {
	  Ioss::Map::build_reverse_map(&reverseNodeMap, &nodeMap[1], nodeCount, 0,
				       "node", myProcessor);
	} else {
	  // Sequential map
	  nodeMap[0] = -1;
	}

      } else {
	unsupported("output nodal id map");
      }
    }
    return nodeMap;
  }

  const Ioss::MapContainer& DatabaseIO::get_element_map() const
  {
    // Allocate space for elemente number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (elementMap.empty()) {
      elementMap.resize(elementCount+1);

      if (is_input()) {
	int error = im_ex_get_elem_num_map(get_file_pointer(), &elementMap[1]);
	if (error < 0) {
	  // Clear out the vector...
	  Ioss::MapContainer().swap(elementMap);
	  pamgen_error(get_file_pointer(), __LINE__, myProcessor);
	}

	// Check for sequential element map.
	// If not, build the reverse G2L element map...
	sequentialEG2L = true;
	for (int i=1; i < elementCount+1; i++) {
	  if (i != elementMap[i]) {
	    sequentialEG2L = false;
	    elementMap[0] = 1;
	    break;
	  }
	}

	if (!sequentialEG2L) {
	  Ioss::Map::build_reverse_map(&reverseElementMap,
				       &elementMap[1], elementCount, 0,
				       "element", myProcessor);
	} else {
	  // Sequential map
	  sequentialEG2L = true;
	  elementMap[0] = -1;
	}

      } else {
	unsupported("output element map");
      }
    }
    return elementMap;
  }

  void DatabaseIO::get_block_adjacencies(const Ioss::ElementBlock *eb,
					 std::vector<std::string> &block_adjacency) const
  {
    if (!blockAdjacenciesCalculated) {
      compute_block_adjacencies();
    }

    // Extract the computed block adjacency information for this
    // element block:
    // Debug print...
    int blk_position  = eb->get_property("original_block_order").get_int();

    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    for (int iblk = 0; iblk < elementBlockCount; iblk++) {

      Ioss::ElementBlock *leb = element_blocks[iblk];
      int lblk_position  = leb->get_property("original_block_order").get_int();

      if (blk_position != lblk_position &&
	  blockAdjacency[blk_position][lblk_position] == 1) {
	block_adjacency.push_back(leb->name());
      }
    }
  }

  void DatabaseIO::compute_block_adjacencies() const
  {
    // Add a field to each element block specifying which other element
    // blocks the block is adjacent to (defined as sharing nodes).
    // This is only calculated on request...

    blockAdjacenciesCalculated = true;
    
    if (elementBlockCount == 1) {
      blockAdjacency.resize(elementBlockCount);
      blockAdjacency[0].resize(elementBlockCount);
      blockAdjacency[0][0] = 0;
      return;
    }

    // Each processor first calculates the adjacencies on their own
    // processor...

    std::vector<int> node_used(nodeCount);
    std::vector<std::vector<int> > inv_con(nodeCount);
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();

    for (int iblk = 0; iblk < elementBlockCount; iblk++) {
      Ioss::ElementBlock *eb = element_blocks[iblk];
      int blk_position =  eb->get_property("original_block_order").get_int();
      int id =            eb->get_property("id").get_int();
      int element_nodes = eb->get_property("topology_node_count").get_int();
      int my_element_count = eb->get_property("entity_count").get_int();
      if (my_element_count > 0) {
	std::vector<int> conn(my_element_count * element_nodes);
	im_ex_get_elem_conn(get_file_pointer(), id, &conn[0]);
	
	for (int i=0; i < my_element_count * element_nodes; i++) {
	  node_used[conn[i]-1] = blk_position+1;
	}

	for (int i=0; i < nodeCount; i++) {
	  if (node_used[i] == blk_position+1) {
	    inv_con[i].push_back(blk_position);
	  }
	}
      }
    }
    
#ifdef HAVE_MPI    
    if (isParallel) {
      // Get contributions from other processors...
      // Get the communication map...
      Ioss::CommSet  *css = get_region()->get_commset("commset_node");
      std::vector<std::pair<int,int> > proc_node;
      {
	std::vector<int> entity_processor;
	css->get_field_data("entity_processor", entity_processor);
	proc_node.resize(entity_processor.size()/2);
	size_t j=0;
	for (size_t i=0; i < entity_processor.size(); j++, i+=2) {
	  proc_node[j].second = entity_processor[i];
	  proc_node[j].first  = entity_processor[i+1];
	}
      }

      // Now sort by increasing processor number.
      std::sort(proc_node.begin() , proc_node.end());
      
      // Pack the data: global_node_id, bits for each block, ...
      // Use 'int' as basic type...
      size_t id_size  = 1;
      size_t word_size = sizeof(int) * 8;
      size_t bits_size = (elementBlockCount+word_size-1)/word_size;
      std::vector<unsigned> send(proc_node.size() * (id_size + bits_size));
      std::vector<unsigned> recv(proc_node.size() * (id_size + bits_size));
      
      std::vector<int> procs(util().parallel_size());
      size_t offset = 0;
      for (size_t i=0; i < proc_node.size(); i++) {
	int glob_id = proc_node[i].second;
	int proc    = proc_node[i].first;
	procs[proc]++;
	send[offset++] = glob_id;
	int loc_id = node_global_to_local(glob_id, true) -1;
	for (size_t j = 0; j < inv_con[loc_id].size(); j++) {
	  int jblk = inv_con[loc_id][j];
	  size_t wrd_off = jblk / word_size;
	  size_t bit     = jblk % word_size;
	  send[offset+wrd_off] |= (1 << bit);
	}
	offset += bits_size;
      }
      
      // Count nonzero entries in 'procs' array -- count of
      // sends/receives
      size_t non_zero = util().parallel_size() - std::count(procs.begin(), procs.end(), 0);
      
      // Post all receives...
      MPI_Request request_null = MPI_REQUEST_NULL;
      std::vector<MPI_Request> request(non_zero, request_null);
      std::vector<MPI_Status>  status(non_zero);
      
      int result = MPI_SUCCESS;
      size_t req_cnt = 0;
      offset = 0;
      for (int i = 0 ; result == MPI_SUCCESS && i < util().parallel_size(); i++) {
	if (procs[i] > 0) {
	  const unsigned size = procs[i] * (id_size + bits_size);
	  void * const   recv_buf  = &recv[offset];
	  result = MPI_Irecv(recv_buf, size, MPI_INT,
			     i, 10101, util().communicator(), &request[req_cnt]);
	  req_cnt++;
	  offset += size;
	}
      }
      assert(result != MPI_SUCCESS || non_zero == req_cnt);

      if (result != MPI_SUCCESS) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Irecv error on processor " << util().parallel_rank()
	       << " in Ioex::DatabaseIO::compute_block_adjacencies";
	std::cerr << errmsg.str();
      }

      int local_error = (MPI_SUCCESS == result) ? 0 : 1 ;
      int global_error = util().global_minmax(local_error, Ioss::ParallelUtils::DO_MAX);

      if (global_error != 0) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Irecv error on some processor "
	       << "in Ioex::DatabaseIO::compute_block_adjacencies";
	IOSS_ERROR(errmsg);
      }
      
      result = MPI_SUCCESS;
      req_cnt = 0;
      offset = 0;
      for (int i = 0 ; result == MPI_SUCCESS && i < util().parallel_size(); i++) {
 	if (procs[i] > 0) {
	  const unsigned size = procs[i] * (id_size + bits_size);
	  void * const   send_buf  = &send[offset];
	  result = MPI_Rsend(send_buf, size, MPI_INT,
			     i, 10101, util().communicator());
	  req_cnt++;
	  offset += size;
	}
      }
      assert(result != MPI_SUCCESS || non_zero == req_cnt);
      
      if (result != MPI_SUCCESS) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Rsend error on processor " << util().parallel_rank()
	       << " in Ioex::DatabaseIO::compute_block_adjacencies";
	std::cerr << errmsg.str();
      }

      local_error = (MPI_SUCCESS == result) ? 0 : 1 ;
      global_error = util().global_minmax(local_error, Ioss::ParallelUtils::DO_MAX);

      if (global_error != 0) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Rsend error on some processor "
	       << "in Ioex::DatabaseIO::compute_block_adjacencies";
	IOSS_ERROR(errmsg);
      }
      
      result = MPI_Waitall(req_cnt, &request[0], &status[0]);
      
      if (result != MPI_SUCCESS) {
	std::ostringstream errmsg;
	errmsg << "ERROR: MPI_Waitall error on processor " << util().parallel_rank()
	       << " in Ioex::DatabaseIO::compute_block_adjacencies";
	std::cerr << errmsg.str();
      }

      // Unpack the data and update the inv_con arrays for boundary
      // nodes...
      offset = 0;
      for (size_t i=0; i < proc_node.size(); i++) {
	int glob_id = recv[offset++];
	int loc_id = node_global_to_local(glob_id, true) -1;
	for (int iblk = 0; iblk < elementBlockCount; iblk++) {
	  size_t wrd_off = iblk / word_size;
	  size_t bit     = iblk % word_size;
	  if (recv[offset+wrd_off] & (1 << bit)) {
	    inv_con[loc_id].push_back(iblk); // May result in duplicates, but that is OK.
	  }
	}
	offset += bits_size;
      }
    }
#endif    
    
    // Convert from inv_con arrays to block adjacency...
    blockAdjacency.resize(elementBlockCount);
    for (int iblk = 0; iblk < elementBlockCount; iblk++) {
      blockAdjacency[iblk].resize(elementBlockCount);
    }

    for (int i = 0; i < nodeCount; i++) {
      for (size_t j = 0; j < inv_con[i].size(); j++) {
	int jblk = inv_con[i][j];
	for (size_t k = j+1; k < inv_con[i].size(); k++) {
	  int kblk = inv_con[i][k];
	  blockAdjacency[jblk][kblk] = 1;
	  blockAdjacency[kblk][jblk] = 1;
	}
      }
    }
	
#ifdef HAVE_MPI
    if (isParallel) {
      // Sync across all processors...
      size_t word_size = sizeof(int) * 8;
      size_t bits_size = (elementBlockCount+word_size-1)/word_size;

      std::vector<unsigned> data(elementBlockCount * bits_size);
      int offset = 0;
      for (int jblk = 0; jblk < elementBlockCount; jblk++) {
	for (int iblk = 0; iblk < elementBlockCount; iblk++) {
	  if (blockAdjacency[jblk][iblk] == 1) {
	    size_t wrd_off = iblk / word_size;
	    size_t bit     = iblk % word_size;
	    data[offset+wrd_off] |= (1 << bit);
	  }
	}
	offset += bits_size;
      }

      std::vector<unsigned> out_data(elementBlockCount * bits_size);
      MPI_Allreduce((void*)&data[0], &out_data[0], static_cast<int>(data.size()),
		    MPI_UNSIGNED, MPI_BOR, util().communicator());

      offset = 0;
      for (int jblk = 0; jblk < elementBlockCount; jblk++) {
	for (int iblk = 0; iblk < elementBlockCount; iblk++) {
	  if (blockAdjacency[jblk][iblk] == 0) {
	    size_t wrd_off = iblk / word_size;
	    size_t bit     = iblk % word_size;
	    if (out_data[offset+wrd_off] & (1 << bit)) {
	      blockAdjacency[jblk][iblk] = 1;
	    }
	  }
	}
	offset += bits_size;
      }
    }
#endif
    
    // Make it symmetric...
    for (int iblk = 0; iblk < elementBlockCount; iblk++) {
      for (int jblk = iblk; jblk < elementBlockCount; jblk++) {
        blockAdjacency[jblk][iblk] = blockAdjacency[iblk][jblk];
      }
    }
  }

  void DatabaseIO::compute_block_membership(int id, std::vector<std::string> &block_membership) const
  {
    Ioss::IntVector block_ids(elementBlockCount);
    if (elementBlockCount == 1) {
      block_ids[0] = 1;
    } else {
      int number_sides;
      int number_distribution_factors;
      int error = im_ex_get_side_set_param(get_file_pointer(), id,
					   &number_sides, &number_distribution_factors);
      if (error < 0) {
	pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }
      
      if (number_sides > 0) {
	// Get the element and element side lists.
	Ioss::IntVector element(number_sides);
	Ioss::IntVector sides(number_sides);
	
	int ierr = im_ex_get_side_set(get_file_pointer(), id, &element[0], &sides[0]);
	if (ierr < 0)
	  pamgen_error(get_file_pointer(), __LINE__, myProcessor);
	
	Ioss::ElementBlock *block = NULL;
	for (int iel = 0; iel < number_sides; iel++) {
	  int elem_id = element[iel];
	  if (block == NULL || !block->contains(elem_id)) {
	    block = get_region()->get_element_block(elem_id);
	    assert(block != NULL);
	    int block_order = block->get_property("original_block_order").get_int();
	    block_ids[block_order] = 1;
	  }
	}
      }

      // Synchronize among all processors....
      if (isParallel) {
	util().global_array_minmax(&block_ids[0],
					 block_ids.size(),
					 Ioss::ParallelUtils::DO_MAX);
      }
    }
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();

    for (int i=0; i < elementBlockCount; i++) {
      if (block_ids[i] == 1) {
	Ioss::ElementBlock *block = element_blocks[i];
	block_membership.push_back(block->name());
      }
    }
  }
  
  void DatabaseIO::compute_block_membership(Ioss::EntityBlock *sideblock,
					    std::vector<std::string> &block_membership) const
  {
    Ioss::IntVector block_ids(elementBlockCount);
    if (elementBlockCount == 1) {
      block_ids[0] = 1;
    } else {
      Ioss::IntVector element_side;
      sideblock->get_field_data("element_side", element_side);
      
      size_t number_sides = element_side.size() / 2;
      Ioss::ElementBlock *block = NULL;
      for (size_t iel = 0; iel < number_sides; iel++) {
	int elem_id = element_side[2*iel];  // Vector contains both element and side.
	elem_id = element_global_to_local(elem_id);
	if (block == NULL || !block->contains(elem_id)) {
	  block = get_region()->get_element_block(elem_id);
	  assert(block != NULL);
	  int block_order = block->get_property("original_block_order").get_int();
	  block_ids[block_order] = 1;
	}
      }
    }
    
    // Synchronize among all processors....
    if (isParallel) {
      util().global_array_minmax(&block_ids[0],
				 block_ids.size(),
				 Ioss::ParallelUtils::DO_MAX);
    }
    
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    
    for (int i=0; i < elementBlockCount; i++) {
      if (block_ids[i] == 1) {
	Ioss::ElementBlock *block = element_blocks[i];
	block_membership.push_back(block->name());
      }
    }
  }

  int DatabaseIO::get_side_connectivity(const Ioss::EntityBlock* fb,
					int id, int,
					int *fconnect,
					size_t /* data_size */) const
  {
    // Get size of data stored on the file...
    int ierr;
    int number_sides;
    int number_distribution_factors;
    ierr = im_ex_get_side_set_param(get_file_pointer(), id,
				    &number_sides, &number_distribution_factors);
    if (ierr < 0)
      pamgen_error(get_file_pointer(), __LINE__, myProcessor);

    // Allocate space for element and local side number
    assert(number_sides > 0);
//----
    Ioss::IntVector element(number_sides);
    Ioss::IntVector side(number_sides);

    ierr = im_ex_get_side_set(get_file_pointer(), id, &element[0], &side[0]);
    if (ierr < 0)
      pamgen_error(get_file_pointer(), __LINE__, myProcessor);
//----

    Ioss::IntVector is_valid_side(number_sides);
    Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, &element[0], &side[0],
						number_sides, get_region());

    Ioss::IntVector elconnect;
    int elconsize = 0; // Size of currently allocated connectivity block
    Ioss::ElementBlock *conn_block = NULL; // Block that we currently
    // have connectivity for

    Ioss::ElementBlock *block = NULL;

    Ioss::MapContainer side_elem_map; // Maps the side into the elements
    // connectivity array
    int current_side = -1;
    int nelnode = 0;
    int nfnodes = 0;
    int ieb = 0;
    int offset = 0;
    for (int iel = 0; iel < number_sides; iel++) {
      if (is_valid_side[iel] == 1) {

	int elem_id = element[iel];

	// ensure we have correct connectivity
	block = get_region()->get_element_block(elem_id);
	if (conn_block != block) {
	  int nelem = block->get_property("entity_count").get_int();
	  nelnode = block->topology()->number_nodes();
	  // Used to map element number into position in connectivity array.
	  // E.g., element 97 is the (97-offset)th element in this block and
	  // is stored in array index (97-offset-1).
	  offset = block->get_offset() + 1;
	  if (elconsize < nelem * nelnode) {
	    elconsize = nelem * nelnode;
	    elconnect.resize(elconsize);
	  }
	  get_field_internal(block, block->get_field("connectivity"),
			     &elconnect[0], nelem*nelnode*sizeof(int));
	  conn_block = block;
	  current_side = -1;
	}

	// NOTE: Element connectivity is returned with nodes in global id space
	if (current_side != side[iel]) {
	  side_elem_map = block->topology()->boundary_connectivity(side[iel]);
	  current_side = side[iel];
	  nfnodes = block->topology()->boundary_type(side[iel])->number_nodes();
	}
	for (int inode = 0; inode < nfnodes; inode++) {
	  int global_node = elconnect[(elem_id-offset)*nelnode +
				      side_elem_map[inode]];
	  fconnect[ieb++] = global_node;
	}
      }
    }
    return ierr;
  }

  // Get distribution factors for the specified side block
  int DatabaseIO::get_side_distributions(const Ioss::EntityBlock* fb,
					 int id, int my_side_count,
					 double *dist_fact,
					 size_t /* data_size */) const
  {
    // NOTE: In pamgen, all distribution factors are 1.0
    const Ioss::ElementTopology *ftopo = fb->topology();
    int nfnodes = ftopo->number_nodes();
    for (int i=0; i < nfnodes * my_side_count; i++)
      dist_fact[i] = 1.0;
    return 0;
  }

  //------------------------------------------------------------------------
}
namespace {
void separate_surface_element_sides(Ioss::IntVector &element,
				    Ioss::IntVector &sides,
				    Ioss::Region *region,
				    Iopg::TopologyMap &topo_map,
				    Iopg::TopologyMap &side_map,
				    Ioss::SurfaceSplitType split_type)
{
  if (element.size() > 0) {
    Ioss::ElementBlock *block = NULL;
    // Topology of sides in current element block
    const Ioss::ElementTopology *common_ftopo = NULL;
    const Ioss::ElementTopology *topo = NULL; // Topology of current side
    int current_side = -1;
    int my_side_count = 0;

    for (size_t iel = 0; iel < element.size(); iel++) {
      int elem_id = element[iel];
      if (block == NULL || !block->contains(elem_id)) {
	block = region->get_element_block(elem_id);
	assert(block != NULL);
	my_side_count = block->topology()->number_boundaries();

	// NULL if hetero sides on element
	common_ftopo = block->topology()->boundary_type(0);
	if (common_ftopo != NULL)
	  topo = common_ftopo;
	current_side = -1;
      }

      if (common_ftopo == NULL && sides[iel] != current_side) {
	current_side = sides[iel];
	assert(current_side > 0 && current_side <= my_side_count);
	topo = block->topology()->boundary_type(sides[iel]);
	assert(topo != NULL);
      }
      std::pair<std::string, const Ioss::ElementTopology*> name_topo;
      if (split_type == Ioss::SPLIT_BY_TOPOLOGIES) {
	name_topo = std::make_pair(block->topology()->name(), topo);
      } else if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
	name_topo = std::make_pair(block->name(), topo);
      }
      topo_map[name_topo]++;
      if (side_map[name_topo] == 0)
	side_map[name_topo] = sides[iel];
      else if (side_map[name_topo] != sides[iel]) {
	// Not a consistent side for all sides in this
	// sideset. Set to large number. Note that maximum
	// sides/element is 6, so don't have to worry about
	// a valid element having 999 sides (unless go to
	// arbitrary polyhedra some time...) Using a large
	// number instead of -1 makes it easier to check the
	// parallel consistency...
	side_map[name_topo] = 999;
      }
    }
  }
}

}
