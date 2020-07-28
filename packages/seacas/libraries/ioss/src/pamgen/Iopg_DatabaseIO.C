// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <create_inline_mesh.h>
#include <pamgen_im_exodusII.h>
#include <pamgen_im_ne_nemesisI.h>

#include <pamgen/Iopg_DatabaseIO.h>

#include <Ioss_CodeTypes.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <float.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits.h>
#include <map>
#include <set>
#include <string>
#include <time.h>
#include <vector>

namespace Iopg {
  using SideSetSet = std::set<std::string>;
  using SideSetMap = std::map<std::string, const std::string, std::less<const std::string>>;

  struct TopologyMapCompare
  {
    bool operator()(const std::pair<std::string, const Ioss::ElementTopology *> &lhs,
                    const std::pair<std::string, const Ioss::ElementTopology *> &rhs) const
    {
      assert(lhs.second != nullptr);
      assert(rhs.second != nullptr);
      return lhs.first < rhs.first ||
             (!(rhs.first < lhs.first) && lhs.second->name() < rhs.second->name());
    }
  };

  using TopologyMap =
      std::map<std::pair<std::string, const Ioss::ElementTopology *>, int, TopologyMapCompare>;
} // namespace Iopg

namespace {
  // Output a message that the operation is unsupported and die...
  void unsupported(const char *operation)
  {
    std::ostringstream errmsg;
    errmsg << "ERROR: Unsupported functionality called: " << operation << '\n';
    IOSS_ERROR(errmsg);
  }

  const size_t max_string_length = MAX_STR_LENGTH;
  const size_t max_line_length   = MAX_LINE_LENGTH;

  void separate_surface_element_sides(Ioss::IntVector &element, Ioss::IntVector &sides,
                                      Ioss::Region *region, Iopg::TopologyMap &topo_map,
                                      Iopg::TopologyMap &    side_map,
                                      Ioss::SurfaceSplitType split_type);

  const char *Version() { return "Iopg_DatabaseIO.C 2010/09/22"; }

  void pamgen_error(int exoid, int lineno, int /* processor */)
  {
    std::ostringstream errmsg;
    errmsg << "Pamgen error at line " << lineno << " in file '" << Version()
           << "' Please report to gdsjaar@sandia.gov if you need help.";
    IOSS_ERROR(errmsg);
  }
} // namespace

namespace Iopg {
  // ========================================================================
  const IOFactory *IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory() : Ioss::IOFactory("pamgen") {}

  Ioss::DatabaseIO *IOFactory::make_IO(const std::string &filename, Ioss::DatabaseUsage db_usage,
                                       MPI_Comm                     communicator,
                                       const Ioss::PropertyManager &properties) const
  {
    return new DatabaseIO(nullptr, filename, db_usage, communicator, properties);
  }

  // ========================================================================
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string &filename,
                         Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
                         const Ioss::PropertyManager &props)
      : Ioss::DatabaseIO(region, filename, db_usage, communicator, props)
  {
    if (is_input()) {
      dbState = Ioss::STATE_UNKNOWN;
    }
    else {
      std::ostringstream errmsg;
      errmsg << "Pamgen mesh option is only valid for input mesh.";
      IOSS_ERROR(errmsg);
    }
  }

  DatabaseIO::~DatabaseIO() { Delete_Pamgen_Mesh(); }

  std::string massagePamgenInputString(std::string sinput, int &dimension)
  {
    std::string        mesh_description;
    std::string        line;
    std::istringstream input(sinput);
    while (std::getline(input, line)) {
      if (line.empty() || line[0] == '$' || line[0] == '#') {
        continue;
      }
      else {
        bool twod = line.find('2') != std::string::npos;
        bool trid = line.find('3') != std::string::npos;
        if (twod && trid) {
          std::ostringstream errmsg;
          errmsg << "PAMGEN: Cannot determine spatial dimension; found both 2 and 3 on first line "
                    "of file.";
          IOSS_ERROR(errmsg);
        }
        else if (twod) {
          dimension = 2;
        }
        else if (trid) {
          dimension = 3;
        }
        break;
      }
    }

    mesh_description += "mesh\n";
    while (std::getline(input, line)) {
      if (line.empty() || line[0] == '$' || line[0] == '#' ||
          (line.find("mesh") != std::string::npos) || (line.find("MESH") != std::string::npos)) {
        ;
      }
      else {
        mesh_description += line;
        mesh_description += "\n";
      }
    }
    return mesh_description;
  }

  void DatabaseIO::read_meta_data__()
  {
    // The file for pamgen contains the mesh description.
    // The Iopg routine is expecting the mesh description to be a
    // single string...

    // Read the data and convert it to a string which is then passed
    // in to the Create_Pamgen_Mesh routine.

    // Check if this is a multi line filename implies it is actually the mesh specification
    std::string mesh_description;

    std::string tfilename = get_filename();
    std::string raw_string;
    // If this is a multi-line string then it actually includes the
    // mesh description; otherwise, it is a file that must be opened
    // and read.
    std::string nlstring = std::string("\n");
    std::size_t found    = tfilename.find(nlstring);
    if (found != std::string::npos) {
      std::istringstream f(tfilename);
      std::string        line;
      while (std::getline(f, line)) {
        raw_string += line;
        raw_string += "\n";
      }
    }
    else {
      std::ifstream f(get_filename());
      if (!f) {
        std::ostringstream errmsg;
        errmsg << "Error opening file '" << get_filename() << "'.";
        IOSS_ERROR(errmsg);
      }
      else {
        std::string line;
        while (std::getline(f, line)) {
          raw_string += line;
          raw_string += "\n";
        }
      }
    }

    // pulls out dimension, gets rid of comments, prepends mesh...
    int dimension    = 3;
    mesh_description = massagePamgenInputString(raw_string, dimension);

    int  retval         = ERROR_FREE_CREATION;
    bool error_detected = false;
    try {
      retval = Create_Pamgen_Mesh(mesh_description.c_str(), dimension, util().parallel_rank(),
                                  util().parallel_size(), INT_MAX);
    }
    catch (const std::exception &x) {
      error_detected = true;
    }

    if (error_detected || retval != ERROR_FREE_CREATION) {
      char *error = nullptr;
      if (retval == ERROR_CREATING_MS) {
        int err_size = getPamgenErrorStreamSize();
        error        = new char[err_size + 1];
        getPamgenErrorStream(error);
      }
      else {
        int err_size = getPamgenEchoStreamSize();
        error        = new char[err_size + 1];
        getPamgenEchoStream(error);
      }
      // handle error
      std::ostringstream errmsg;
      errmsg << "Pamgen mesh generation failed:\n";
      errmsg << error;
      delete[] error;
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

    char dbtitle[max_line_length + 1];
    std::memset(dbtitle, 0, max_line_length + 1);

    int node_count = 0;
    int elem_count = 0;
    int error      = im_ex_get_init(get_file_pointer(), dbtitle, &spatialDimension, &node_count,
                               &elem_count, &elementBlockCount, &nodesetCount, &sidesetCount);
    if (error < 0)
      pamgen_error(get_file_pointer(), __LINE__, myProcessor);

    nodeCount    = node_count;
    elementCount = elem_count;

    nodeBlockCount = 1;

    if (nodeCount == 0) {
      Ioss::WARNING() << "No nodes were found in the model, file '" << decoded_filename() << "'";
    }
    else if (nodeCount < 0) {
      // NOTE: Code will not continue past this call...
      std::ostringstream errmsg;
      errmsg << "Negative node count was found in the model\n"
             << "       File: '" << decoded_filename() << "'.\n";
      IOSS_ERROR(errmsg);
    }

    if (elementCount == 0) {
      Ioss::WARNING() << "No elements were found in the model, file: '" << decoded_filename()
                      << "'";
    }

    if (elementCount < 0) {
      // NOTE: Code will not continue past this call...
      std::ostringstream errmsg;
      errmsg << "Negative element count was found in the model, file: '" << decoded_filename()
             << "'";
      IOSS_ERROR(errmsg);
    }

    if (elementCount > 0 && elementBlockCount <= 0) {
      // NOTE: Code will not continue past this call...
      std::ostringstream errmsg;
      errmsg << "No element blocks were found in the model, file: '" << decoded_filename() << "'";
      IOSS_ERROR(errmsg);
    }

    Ioss::Region *this_region = get_region();
    this_region->property_add(Ioss::Property(std::string("title"), dbtitle));
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
    int  num_proc;         // Number of processors file was decomposed for
    int  num_proc_in_file; // Number of processors this file has info for
    char file_type[2];     // "s" for scalar, "p" for parallel
    // Sierra is set up to handle 1 processor per file; parallel file

    int error =
        im_ne_get_init_info(get_file_pointer(), &num_proc, &num_proc_in_file, &file_type[0]);
    if (error < 0) {
      // Not a nemesis file
      if (util().parallel_size() > 1) {
        std::ostringstream errmsg;
        errmsg << "Exodus file does not contain nemesis information.\n";
        IOSS_ERROR(errmsg);
      }
    }
    else if (num_proc != util().parallel_size() && util().parallel_size() > 1) {
      std::ostringstream errmsg;
      errmsg << "Exodus file was decomposed for " << num_proc
             << " processors; application is currently being run on " << util().parallel_size()
             << " processors";
      IOSS_ERROR(errmsg);
    }
    else if (num_proc_in_file != 1) {
      std::ostringstream errmsg;
      errmsg << "Exodus file contains data for " << num_proc_in_file
             << " processors; application requires 1 processor per file.";
      IOSS_ERROR(errmsg);
    }
    else if (file_type[0] != 'p') {
      std::ostringstream errmsg;
      errmsg << "Exodus file contains scalar nemesis data; application requires parallel nemesis "
                "data.";
      IOSS_ERROR(errmsg);
    }

    // Get global data (over all processors)
    int global_nodes    = nodeCount;
    int global_elements = elementCount;
    int global_eblocks  = 0; // unused
    int global_nsets    = 0; // unused
    int global_ssets    = 0; // unused

    int num_external_nodes; // unused
    int num_elem_cmaps     = 0;
    int num_node_cmaps     = 0;
    int num_internal_nodes = nodeCount;
    int num_border_nodes   = 0;
    int num_internal_elems = elementCount;
    int num_border_elems   = 0;

    if (isParallel) {

      error = im_ne_get_loadbal_param(get_file_pointer(), &num_internal_nodes, &num_border_nodes,
                                      &num_external_nodes, &num_internal_elems, &num_border_elems,
                                      &num_node_cmaps, &num_elem_cmaps, myProcessor);
      if (error < 0)
        pamgen_error(get_file_pointer(), __LINE__, myProcessor);

      commsetNodeCount = num_node_cmaps;
      commsetElemCount = num_elem_cmaps;

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
    region->property_add(Ioss::Property("internal_node_count", num_internal_nodes));
    region->property_add(Ioss::Property("border_node_count", num_border_nodes));
    region->property_add(Ioss::Property("internal_element_count", num_internal_elems));
    region->property_add(Ioss::Property("border_element_count", num_border_elems));
    region->property_add(Ioss::Property("global_node_count", global_nodes));
    region->property_add(Ioss::Property("global_element_count", global_elements));

    // Possibly, the following 4 fields should be nodesets and element
    // sets instead of fields on the region...
    region->field_add(Ioss::Field("internal_nodes", Ioss::Field::INTEGER, IOSS_SCALAR(),
                                  Ioss::Field::COMMUNICATION, num_internal_nodes));
    region->field_add(Ioss::Field("border_nodes", Ioss::Field::INTEGER, IOSS_SCALAR(),
                                  Ioss::Field::COMMUNICATION, num_border_nodes));
    region->field_add(Ioss::Field("internal_elements", Ioss::Field::INTEGER, IOSS_SCALAR(),
                                  Ioss::Field::COMMUNICATION, num_internal_elems));
    region->field_add(Ioss::Field("border_elements", Ioss::Field::INTEGER, IOSS_SCALAR(),
                                  Ioss::Field::COMMUNICATION, num_border_elems));

    assert(nodeCount == num_internal_nodes + num_border_nodes);
    assert(elementCount == num_internal_elems + num_border_elems);
  }

  void DatabaseIO::get_nodeblocks()
  {
    // For exodusII, there is only a single node block which contains
    // all of the nodes.
    // The default id assigned is '1' and the name is 'nodeblock_1'
    std::string      block_name = "nodeblock_1";
    Ioss::NodeBlock *block = new Ioss::NodeBlock(this, block_name, nodeCount, spatialDimension);
    block->property_add(Ioss::Property("id", 1));
    block->property_add(Ioss::Property("guid", util().generate_guid(1)));

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

    size_t            all_element_type_length = elementBlockCount * (max_string_length + 1);
    std::vector<char> all_element_type(all_element_type_length);

    Ioss::IntVector attributes(elementBlockCount);
    Ioss::IntVector my_node_count(elementBlockCount);
    Ioss::IntVector local_element_count(elementBlockCount);
    Ioss::IntVector global_element_count(elementBlockCount);
    int             iblk;

    for (iblk = 0; iblk < elementBlockCount; iblk++) {
      int id = element_block_ids[iblk];
      int number_elements;
      int nodes_per_element;
      int attributes_per_element;

      char *const element_type = &all_element_type[0] + iblk * (max_string_length + 1);

      error = im_ex_get_elem_block(get_file_pointer(), id, element_type, &number_elements,
                                   &nodes_per_element, &attributes_per_element);
      if (error < 0) {
        pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }

      local_element_count[iblk] = number_elements;

      if (number_elements == 0) {
        attributes[iblk]    = 0;
        my_node_count[iblk] = 0;
      }
      else {
        attributes[iblk]    = attributes_per_element;
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

    for (iblk = 0; iblk < elementBlockCount; iblk++) {
      int         id           = element_block_ids[iblk];
      std::string alias        = Ioss::Utils::encode_entity_name("block", id);
      char *const element_type = &all_element_type[0] + iblk * (max_string_length + 1);

      Ioss::ElementBlock *block      = nullptr;
      std::string         block_name = Ioss::Utils::encode_entity_name("block", id);

      Ioss::Utils::fixup_name(element_type); // Convert to lowercase; replace spaces with '_'
      std::string type      = std::string(element_type);
      std::string save_type = type;

      // Fixup an exodusII kluge/ambiguity.
      // The element block type does not fully define the element. For
      // example, a block of type 'triangle' may have either 3 or 6
      // nodes.  To fix this, check the block type name and see if it
      // ends with a number.  If it does, assume it is OK; if not, append
      // the 'nodes_per_element'.
      if (!isdigit(*(type.rbegin()))) {
        if (my_node_count[iblk] > 1) {
          type += std::to_string(my_node_count[iblk]);
        }
      }

      block = new Ioss::ElementBlock(this, block_name, type, local_element_count[iblk]);

      block->property_add(Ioss::Property("id", id));
      block->property_add(Ioss::Property("guid", util().generate_guid(id)));
      block->property_add(Ioss::Property("original_block_order", iblk));

      if (block->get_property("topology_type").get_string() != save_type && save_type != "null" &&
          save_type != "") {
        // Maintain original element type on output database if possible.
        block->property_add(Ioss::Property("original_topology_type", save_type));
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
      int             error = im_ex_get_node_set_ids(get_file_pointer(), &nodeset_ids[0]);
      if (error < 0) {
        pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }

      for (int ins = 0; ins < nodesetCount; ins++) {
        int id = nodeset_ids[ins];
        int number_nodes;
        int number_distribution_factors;

        error = im_ex_get_node_set_param(get_file_pointer(), id, &number_nodes,
                                         &number_distribution_factors);
        if (error < 0) {
          pamgen_error(get_file_pointer(), __LINE__, myProcessor);
        }

        std::string    nodeset_name = Ioss::Utils::encode_entity_name("nodelist", id);
        Ioss::NodeSet *nodeset      = new Ioss::NodeSet(this, nodeset_name, number_nodes);
        nodeset->property_add(Ioss::Property("id", id));
        nodeset->property_add(Ioss::Property("guid", util().generate_guid(id)));
        get_region()->add(nodeset);

        get_region()->add_alias(nodeset_name, Ioss::Utils::encode_entity_name("nodelist", id));
        get_region()->add_alias(nodeset_name, Ioss::Utils::encode_entity_name("nodeset", id));
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
      int elem_count    = 0;

      // NOTE: It is possible for a parallel run to have no
      // communications maps if the decomposition occurs along contact
      // surfaces.  In this case, we create empty node and element
      // communication maps.

      if (commsetNodeCount > 0 || commsetElemCount > 0) {
        if (commsetNodeCount > 0) {
          nodeCmapIds.resize(commsetNodeCount);
          nodeCmapNodeCnts.resize(commsetNodeCount);
        }
        if (commsetElemCount > 0) {
          elemCmapIds.resize(commsetElemCount);
          elemCmapElemCnts.resize(commsetElemCount);
        }

        int error =
            im_ne_get_cmap_params(get_file_pointer(), nodeCmapIds.data(), nodeCmapNodeCnts.data(),
                                  elemCmapIds.data(), elemCmapElemCnts.data(), myProcessor);
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
      Ioss::CommSet *commset = new Ioss::CommSet(this, "commset_node", "node", my_node_count);
      commset->property_add(Ioss::Property("id", 1));
      commset->property_add(Ioss::Property("guid", util().generate_guid(1)));
      get_region()->add(commset);

      commset = new Ioss::CommSet(this, "commset_side", "side", elem_count);
      commset->property_add(Ioss::Property("id", 1));
      commset->property_add(Ioss::Property("guid", util().generate_guid(1)));
      get_region()->add(commset);
    }
  }

  void DatabaseIO::get_sidesets()
  {
    // This function creates all sidesets (surfaces) for a
    // model.  Note that a side set contains 1 or more side
    // blocks which are homogeneous (same topology). In serial execution,
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
      int             error = im_ex_get_side_set_ids(get_file_pointer(), &side_set_ids[0]);
      if (error < 0) {
        pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }

      for (int iss = 0; iss < sidesetCount; iss++) {
        int         id  = side_set_ids[iss];
        std::string sid = "";
        int         number_sides;
        int         number_distribution_factors;
        TopologyMap topo_map;
        TopologyMap side_map; // Used to determine side consistency

        Ioss::SurfaceSplitType split_type = splitType;
        std::string            side_set_name;

        side_set_name = Ioss::Utils::encode_entity_name("surface", id);

        Ioss::SideSet *side_set = new Ioss::SideSet(this, side_set_name);
        get_region()->add(side_set);
        side_set->property_add(Ioss::Property("id", id));
        side_set->property_add(Ioss::Property("guid", util().generate_guid(id)));

        get_region()->add_alias(side_set_name, Ioss::Utils::encode_entity_name("surface", id));
        get_region()->add_alias(side_set_name, Ioss::Utils::encode_entity_name("sideset", id));

        //        split_type = SPLIT_BY_ELEMENT_BLOCK;
        //        split_type = SPLIT_BY_TOPOLOGIES;
        //        split_type = SPLIT_BY_DONT_SPLIT;

        // Determine how many side blocks compose this side set.
        error = im_ex_get_side_set_param(get_file_pointer(), id, &number_sides,
                                         &number_distribution_factors);
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
          topo_map[std::make_pair(sideTopology[0].first->name(), sideTopology[0].second)] =
              number_sides;
        }
        else if (split_type == Ioss::SPLIT_BY_DONT_SPLIT) {
          const Ioss::ElementTopology *mixed_topo = Ioss::ElementTopology::factory("unknown");
          topo_map[std::make_pair("unknown", mixed_topo)] = number_sides;
        }
        else if (split_type == Ioss::SPLIT_BY_TOPOLOGIES) {
          // There are multiple side types in the model.
          // Iterate through the elements in the sideset, determine
          // their parent element block using the blocks element
          // topology and the side number, determine the side
          // type.

          for (Ioss::TopoContainer::size_type i = 0; i < sideTopology.size(); i++) {
            topo_map[std::make_pair(sideTopology[i].first->name(), sideTopology[i].second)] = 0;
            side_map[std::make_pair(sideTopology[i].first->name(), sideTopology[i].second)] = 0;
          }

          separate_surface_element_sides(element, sides, get_region(), topo_map, side_map,
                                         split_type);
        }
        else if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
          // There are multiple side types in the model.  Iterate
          // through the elements in the sideset, determine their parent
          // element block using blocks element topology and the side
          // number, determine the side type.

          // Seed the topo_map map with <block->name, face_topo>
          // pairs so we are sure that all processors have the same
          // starting topo_map (size and order).
          const Ioss::ElementBlockContainer &element_blocks = get_region()->get_element_blocks();

          for (int i = 0; i < elementBlockCount; i++) {
            Ioss::ElementBlock *         block        = element_blocks[i];
            const std::string &          name         = block->name();
            const Ioss::ElementTopology *common_ftopo = block->topology()->boundary_type(0);
            if (common_ftopo != nullptr) {
              // All sides of this element block's topology have the same topology
              topo_map[std::make_pair(name, common_ftopo)] = 0;
              side_map[std::make_pair(name, common_ftopo)] = 0;
            }
            else {
              // The sides have different topology, iterate over
              // them and create an entry for the unique face
              // topology types
              int par_dim = block->topology()->parametric_dimension();
              if (par_dim == 2 || par_dim == 3) {
                int my_side_count = block->topology()->number_boundaries();
                for (int ii = 0; ii < my_side_count; ii++) {
                  const Ioss::ElementTopology *topo    = block->topology()->boundary_type(ii + 1);
                  topo_map[std::make_pair(name, topo)] = 0;
                  side_map[std::make_pair(name, topo)] = 0;
                }
              }
            }
          }
          separate_surface_element_sides(element, sides, get_region(), topo_map, side_map,
                                         split_type);
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
          util().global_array_minmax(global_side_counts, Ioss::ParallelUtils::DO_MAX);
        }

        // Create Side Blocks

        int                         i = 0;
        TopologyMap::const_iterator I = topo_map.begin();
        assert(I != topo_map.end());
        while (I != topo_map.end()) {
          if (global_side_counts[i++] > 0) {
            const std::string            topo_or_block_name = (*I).first.first;
            const Ioss::ElementTopology *side_topo          = (*I).first.second;
            assert(side_topo != nullptr);
#if 0
            if (side_topo->parametric_dimension() == topology_dimension-1 ||
                split_type == Ioss::SPLIT_BY_DONT_SPLIT ) {
#else
            if (true) {
#endif
            int my_side_count = (*I).second;

            std::string side_block_name = "surface_" + topo_or_block_name + "_" + side_topo->name();
            if (side_set_name == "universal_sideset") {
              side_block_name = side_set_name;
            }
            else {
              if (sid == "")
                side_block_name = Ioss::Utils::encode_entity_name(side_block_name, id);
              else {
                side_block_name += "_";
                side_block_name += sid;
              }
            }

            Ioss::ElementBlock *block = nullptr;
            // Need to get elem_topo....
            const Ioss::ElementTopology *elem_topo = nullptr;
            if (split_type == Ioss::SPLIT_BY_TOPOLOGIES) {
              elem_topo = Ioss::ElementTopology::factory(topo_or_block_name);
            }
            else if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
              block = get_region()->get_element_block(topo_or_block_name);
              if (block == nullptr) {
                std::ostringstream errmsg;
                errmsg << "INTERNAL ERROR: Could not find element block '" << topo_or_block_name
                       << "' Something is wrong in the Iopg::DatabaseIO class. Please report.\n";
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
            assert(elem_topo != nullptr);

            Ioss::SideBlock *side_block = new Ioss::SideBlock(
                this, side_block_name, side_topo->name(), elem_topo->name(), my_side_count);
            side_set->add(side_block);

            // Note that all sideblocks within a specific
            // sideset might have the same id.
            assert(side_block != nullptr);
            side_block->property_add(Ioss::Property("id", id));
            side_block->property_add(Ioss::Property("guid", util().generate_guid(id)));

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
            if (!side_map.empty()) {
              // Set a property indicating which element side
              // (1-based) all faces in this block are applied to.
              // If they are not all assigned to the same element
              // side, indicate this with a side equal to 0.
              //
              // (note: 'i' has already been incremented earlier in
              // the loop.  We need previous value here...)
              int side = global_side_counts[i - 1 + topo_map.size()];
              if (side == 999)
                side = 0;
              assert(side <= elem_topo->number_boundaries());
              side_block->set_consistent_side_number(side);
            }

            // Add an alias...
            get_region()->add_alias(side_block);

            if (split_type != Ioss::SPLIT_BY_DONT_SPLIT && side_set_name != "universal_sideset") {
              std::string storage = "Real[";
              storage += std::to_string(side_topo->number_nodes());
              storage += "]";
              side_block->field_add(Ioss::Field("distribution_factors", Ioss::Field::REAL, storage,
                                                Ioss::Field::MESH, my_side_count));
            }

            if (side_set_name == "universal_sideset") {
              side_block->field_add(Ioss::Field("side_ids", Ioss::Field::INTEGER, "scalar",
                                                Ioss::Field::MESH, my_side_count));
            }
          }
        }
        ++I;
      }
    }
  }
} // namespace Iopg

bool DatabaseIO::begin__(Ioss::State /* state */) { return true; }

bool DatabaseIO::end__(Ioss::State /* state */) { return true; }

int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                       void *data, size_t data_size) const
{
  size_t num_to_get = field.verify(data_size);
  if (num_to_get > 0) {

    Ioss::Field::RoleType role = field.get_role();
    if (role == Ioss::Field::MESH) {

      if ((field.get_name() == "mesh_model_coordinates_x") ||
          (field.get_name() == "mesh_model_coordinates_y") ||
          (field.get_name() == "mesh_model_coordinates_z")) {

        std::vector<double> x(num_to_get);
        std::vector<double> y(num_to_get);
        std::vector<double> z;
        if (spatialDimension == 3)
          z.resize(num_to_get);

        double *rdata = static_cast<double *>(data);

        int ierr = im_ex_get_coord(get_file_pointer(), &x[0], &y[0], &z[0]);
        if (ierr < 0) {
          pamgen_error(get_file_pointer(), __LINE__, myProcessor);
        }
        size_t index = 0;

        if (field.get_name() == "mesh_model_coordinates_x") {
          for (size_t i = 0; i < num_to_get; i++) {
            rdata[index++] = x[i];
          }
        }
        else if (field.get_name() == "mesh_model_coordinates_y") {
          for (size_t i = 0; i < num_to_get; i++) {
            rdata[index++] = y[i];
          }
        }
        else {
          for (size_t i = 0; i < num_to_get; i++) {
            if (spatialDimension == 3)
              rdata[index++] = z[i];
          }
        }
      }

      else if (field.get_name() == "mesh_model_coordinates") {
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
        double *rdata = static_cast<double *>(data);

        int ierr = im_ex_get_coord(get_file_pointer(), &x[0], &y[0], &z[0]);
        if (ierr < 0)
          pamgen_error(get_file_pointer(), __LINE__, myProcessor);

        size_t index = 0;
        for (size_t i = 0; i < num_to_get; i++) {
          rdata[index++] = x[i];
          rdata[index++] = y[i];
          if (spatialDimension == 3)
            rdata[index++] = z[i];
        }
      }

      // NOTE: The implicit_ids field is ONLY provided for backward-
      // compatibility and should not be used unless absolutely
      // required. For generated mesh, the implicit_ids and ids are the same.
      else if (field.get_name() == "ids" || field.get_name() == "implicit_ids") {
        // Map the local ids in this node block
        // (1...node_count) to global node ids.
        get_node_map().map_implicit_data(data, field, num_to_get, 0);
      }
      else if (field.get_name() == "connectivity") {
        // Do nothing, just handles an idiosyncrasy of the GroupingEntity
      }
      else if (field.get_name() == "owning_processor") {
        if (isParallel) {
          Ioss::CommSet *css   = get_region()->get_commset("commset_node");
          int *          idata = static_cast<int *>(data);
          for (size_t i = 0; i < num_to_get; i++) {
            idata[i] = myProcessor;
          }

          std::vector<int> ent_proc;
          css->get_field_data("entity_processor_raw", ent_proc);
          for (size_t i = 0; i < ent_proc.size(); i += 2) {
            int node = ent_proc[i + 0];
            int proc = ent_proc[i + 1];
            if (proc < myProcessor) {
              idata[node - 1] = proc;
            }
          }
        }
        else {
          // Serial case...
          int *idata = static_cast<int *>(data);
          for (int64_t i = 0; i < nodeCount; i++) {
            idata[i] = 0;
          }
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(nb, field, "input");
      }
    }
  }
  return num_to_get;
}

int64_t DatabaseIO::get_field_internal(const Ioss::Region * /* reg */,
                                       const Ioss::Field & /* field */, void * /* data */,
                                       size_t /* data_size */) const
{
  return -1;
}

int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                       void *data, size_t data_size) const
{
  size_t num_to_get = field.verify(data_size);
  if (num_to_get > 0) {

    int                   ierr             = 0;
    int                   id               = eb->get_property("id").get_int();
    int                   my_element_count = eb->entity_count();
    Ioss::Field::RoleType role             = field.get_role();

    if (role == Ioss::Field::MESH) {
      // Handle the MESH fields required for an ExodusII file model.
      // (The 'genesis' portion)

      if (field.get_name() == "connectivity" || field.get_name() == "connectivity_raw") {
        int element_nodes = eb->get_property("topology_node_count").get_int();
        assert(field.raw_storage()->component_count() == element_nodes);

        // The connectivity is stored in a 1D array.
        // The element_node index varies fastet
        if (my_element_count > 0) {
          ierr = im_ex_get_elem_conn(get_file_pointer(), id, static_cast<int *>(data));
          if (ierr < 0)
            pamgen_error(get_file_pointer(), __LINE__, myProcessor);

          // Now, map the nodes in the connectivity from local to global ids
          if (field.get_name() == "connectivity") {
            get_node_map().map_data(data, field, num_to_get * element_nodes);
          }
        }
      }
      else if (field.get_name() == "ids") {
        // Map the local ids in this element block
        // (eb_offset+1...eb_offset+1+my_element_count) to global element ids.
        get_element_map().map_implicit_data(data, field, num_to_get, eb->get_offset());
      }
      else {
        num_to_get = Ioss::Utils::field_warning(eb, field, "input");
      }
    }

    else if (role == Ioss::Field::ATTRIBUTE) {
      // ERROR: Pamgen doesn't support elements with attributes...
    }
    else if (role == Ioss::Field::TRANSIENT) {
      // ERROR: Pamgen doesn't support transient data...
    }
    else if (role == Ioss::Field::REDUCTION) {
      num_to_get = Ioss::Utils::field_warning(eb, field, "input reduction");
    }
  }
  return num_to_get;
}

int64_t DatabaseIO::get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                       void *data, size_t data_size) const
{
  {
    Ioss::SerializeIO serializeIO__(this);

    size_t num_to_get = field.verify(data_size);

    if (num_to_get > 0) {
      int entity_count = cs->entity_count();

      // Return the <entity (node or side), processor> pair
      if (field.get_name() == "entity_processor" || field.get_name() == "entity_processor_raw") {

        // Check type -- node or side
        std::string type = cs->get_property("entity_type").get_string();

        // Allocate temporary storage space
        Ioss::IntVector entities(num_to_get);
        Ioss::IntVector procs(num_to_get);

        if (type == "node") {
          int cm_offset = 0;
          int i;
          for (i = 0; i < commsetNodeCount; i++) {
            int ierr = im_ne_get_node_cmap(get_file_pointer(), nodeCmapIds[i], &entities[cm_offset],
                                           &procs[cm_offset], myProcessor);
            if (ierr < 0)
              pamgen_error(get_file_pointer(), __LINE__, myProcessor);
            cm_offset += nodeCmapNodeCnts[i];
          }
          assert(cm_offset == entity_count);

          // Convert local node id to global node id and store in 'data'
          int *entity_proc = static_cast<int *>(data);
          if (field.get_name() == "entity_processor") {
            const Ioss::MapContainer &map = get_node_map().map();

            int j = 0;
            for (i = 0; i < entity_count; i++) {
              int local_id     = entities[i];
              entity_proc[j++] = map[local_id];
              entity_proc[j++] = procs[i];
            }
          }
          else {
            int j = 0;
            for (i = 0; i < entity_count; i++) {
              entity_proc[j++] = entities[i];
              entity_proc[j++] = procs[i];
            }
          }
        }
        else if (type == "side") {
          Ioss::IntVector sides(entity_count);
          int             cm_offset = 0;
          int             i;
          for (i = 0; i < commsetElemCount; i++) {
            int ierr = im_ne_get_elem_cmap(get_file_pointer(), elemCmapIds[i], &entities[cm_offset],
                                           &sides[cm_offset], &procs[cm_offset], myProcessor);
            if (ierr < 0)
              pamgen_error(get_file_pointer(), __LINE__, myProcessor);
            cm_offset += elemCmapElemCnts[i];
          }
          assert(cm_offset == entity_count);

          int *entity_proc = static_cast<int *>(data);
          if (field.get_name() == "entity_processor") {
            const Ioss::MapContainer &map = get_element_map().map();

            int j = 0;
            for (i = 0; i < entity_count; i++) {
              int local_id     = entities[i];
              entity_proc[j++] = map[local_id];
              entity_proc[j++] = sides[i];
              entity_proc[j++] = procs[i];
            }
          }
          else {
            int j = 0;
            for (i = 0; i < entity_count; i++) {
              entity_proc[j++] = entities[i];
              entity_proc[j++] = sides[i];
              entity_proc[j++] = procs[i];
            }
          }
        }
        else {
          std::ostringstream errmsg;
          errmsg << "Invalid commset type " << type;
          IOSS_ERROR(errmsg);
        }
      }
      else if (field.get_name() == "ids") {
        // Do nothing, just handles an idiosyncrasy of the GroupingEntity
      }
      else {
        num_to_get = Ioss::Utils::field_warning(cs, field, "input");
      }
    }
    return num_to_get;
  }
}

int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock *fb, const Ioss::Field &field,
                                       void *data, size_t data_size) const
{
  size_t num_to_get = field.verify(data_size);
  if (num_to_get > 0) {

    int    id           = fb->get_property("id").get_int();
    size_t entity_count = fb->entity_count();
    if (num_to_get != entity_count) {
      std::ostringstream errmsg;
      errmsg << "Partial field input not yet implemented for side blocks";
      IOSS_ERROR(errmsg);
    }

    int number_sides;
    int number_distribution_factors;
    int ierr = im_ex_get_side_set_param(get_file_pointer(), id, &number_sides,
                                        &number_distribution_factors);
    if (ierr < 0)
      pamgen_error(get_file_pointer(), __LINE__, myProcessor);

    Ioss::Field::RoleType role = field.get_role();
    if (role == Ioss::Field::MESH) {

      // In exodusII, we may have split the sideset into multiple
      // side blocks if there are multiple side  topologies in the
      // sideset.  Because of this, the passed in 'data' may not be
      // large enough to hold the data residing in the sideset and we
      // may need to allocate a temporary array...  This can be checked
      // by comparing the size of the sideset with the 'my_side_count' of
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
        const Ioss::MapContainer &map = get_element_map().map();

        // Allocate space for local side number, use 'data' as temporary
        // storage for element numbers and overwrite with the side
        // numbers.
        Ioss::IntVector sides;
        int *           element = nullptr;
        int *           ids     = static_cast<int *>(data);
        if (number_sides == static_cast<int>(entity_count)) {
          // Only 1 side block in this sideset
          sides.resize(entity_count);
          element = static_cast<int *>(data);
        }
        else {
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
            int new_id = 10 * map[element[iel]] + sides[iel];
            ids[iel]   = new_id;
          }
        }
        else {
          Ioss::IntVector is_valid_side;
          Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, 4, element, sides.data(),
                                                      number_sides, get_region());
          size_t ieb = 0;
          for (int iel = 0; iel < number_sides; iel++) {
            if (is_valid_side[iel] == 1) {
              int new_id = 10 * map[element[iel]] + sides[iel];
              ids[ieb++] = new_id;
            }
          }
          assert(ieb == entity_count);
        }

        if (number_sides != static_cast<int>(entity_count))
          delete[] element;
      }
      else if (field.get_name() == "element_side") {
        // In exodusII, the 'side set' is stored as a sideset.  A sideset
        // has a list of elements and a corresponding local element side
        // (1-based)

        // Since we only have a single array, we need to allocate an extra
        // array to store all of the data.  Note also that the element_id
        // is the global id but only the local id is stored so we need to
        // map from local_to_global prior to generating the side  id...

        // Get the element number map (1-based)...
        const Ioss::MapContainer &map = get_element_map().map();

        // Allocate space for local side number and element numbers
        // numbers.
        int *element_side = static_cast<int *>(data);

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
          assert(index / 2 == entity_count);
        }
        else {
          Ioss::IntVector is_valid_side;
          Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, 4, element.data(),
                                                      sides.data(), number_sides, get_region());

          size_t index = 0;
          for (int iel = 0; iel < number_sides; iel++) {
            if (is_valid_side[iel] == 1) {
              // This side  belongs in the side block
              element_side[index++] = map[element[iel]];
              element_side[index++] = sides[iel];
            }
          }
          assert(index / 2 == entity_count);
        }
      }
      else if (field.get_name() == "connectivity") {
        // The side  connectivity needs to be generated 'on-the-fly' from
        // the element number and local side  of that element. A sideset
        // can span multiple element blocks, and contain multiple side
        // types; the side block contains side of similar topology.
        ierr = get_side_connectivity(fb, id, entity_count, static_cast<int *>(data),
                                     data_size / sizeof(int));
        if (ierr < 0)
          pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }
      else if (field.get_name() == "distribution_factors") {
        ierr = get_side_distributions(fb, id, entity_count, static_cast<double *>(data),
                                      data_size / sizeof(double));
        if (ierr < 0)
          pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }
      else {
        num_to_get = Ioss::Utils::field_warning(fb, field, "input");
      }
    }
    else if (role == Ioss::Field::TRANSIENT) {
      unsupported("side set transient fields");
    }
  }
  return num_to_get;
}

int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                       void *data, size_t data_size) const
{
  size_t num_to_get = field.verify(data_size);
  if (num_to_get > 0) {

    int                   id   = ns->get_property("id").get_int();
    Ioss::Field::RoleType role = field.get_role();
    if (role == Ioss::Field::MESH) {

      if (field.get_name() == "ids") {
        int ierr = im_ex_get_node_set(get_file_pointer(), id, static_cast<int *>(data));
        if (ierr < 0)
          pamgen_error(get_file_pointer(), __LINE__, myProcessor);

        // Convert the local node ids to global ids
        get_node_map().map_data(data, field, num_to_get);
      }
      else if (field.get_name() == "distribution_factors") {
        // Not supported by pamgen.  Fill data with '1'
        double *rdata = static_cast<double *>(data);
        for (size_t i = 0; i < num_to_get; i++)
          rdata[i] = 1.0;
      }
      else {
        num_to_get = Ioss::Utils::field_warning(ns, field, "input");
      }
    }
    else if (role == Ioss::Field::TRANSIENT) {
      unsupported("transient nodeset fields");
    }
  }
  return num_to_get;
}

int64_t DatabaseIO::get_field_internal(const Ioss::SideSet * /* fs */,
                                       const Ioss::Field & /* field */, void * /* data */,
                                       size_t /* data_size */) const
{
  return -1;
}
int64_t DatabaseIO::put_field_internal(const Ioss::Region * /* region */,
                                       const Ioss::Field & /* field */, void * /* data */,
                                       size_t /* data_size */) const
{
  unsupported("output region field");
  return -1;
}

int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock * /* eb */,
                                       const Ioss::Field & /* field */, void * /* data */,
                                       size_t /* data_size */) const
{
  unsupported("output element block field");
  return -1;
}
int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock * /* fb */,
                                       const Ioss::Field & /* field */, void * /* data */,
                                       size_t /* data_size */) const
{
  unsupported("output sideblock field");
  return -1;
}
int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock * /* nb */,
                                       const Ioss::Field & /* field */, void * /* data */,
                                       size_t /* data_size */) const
{
  unsupported("output nodeblock field");
  return -1;
}

int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet * /* ns */,
                                       const Ioss::Field & /* field */, void * /* data */,
                                       size_t /* data_size */) const
{
  unsupported("output nodeset field");
  return -1;
}
int64_t DatabaseIO::put_field_internal(const Ioss::SideSet * /* fs */,
                                       const Ioss::Field & /* field */, void * /* data */,
                                       size_t /* data_size */) const
{
  unsupported("output sideset field");
  return -1;
}
int64_t DatabaseIO::put_field_internal(const Ioss::CommSet * /* cs */,
                                       const Ioss::Field & /* field */, void * /* data */,
                                       size_t /* data_size */) const
{
  unsupported("output commset field");
  return -1;
}

const Ioss::Map &DatabaseIO::get_node_map() const
{
  // Allocate space for node number map and read it in...
  // Can be called multiple times, allocate 1 time only
  if (nodeMap.map().empty()) {
    nodeMap.set_size(nodeCount);

    if (is_input()) {
      std::vector<int> node_map(nodeMap.map().size() - 1);
      int              error = im_ex_get_node_num_map(get_file_pointer(), node_map.data());
      if (error < 0) {
        // Clear out the vector...
        Ioss::MapContainer().swap(nodeMap.map());
        pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }
      nodeMap.set_map(node_map.data(), node_map.size(), 0, true);
    }
    else {
      unsupported("output nodal id map");
    }
  }
  return nodeMap;
}

const Ioss::Map &DatabaseIO::get_element_map() const
{
  // Allocate space for element number map and read it in...
  // Can be called multiple times, allocate 1 time only
  if (elemMap.map().empty()) {
    elemMap.set_size(elementCount);

    if (is_input()) {
      std::vector<int> elem_map(elemMap.map().size() - 1);
      int              error = im_ex_get_elem_num_map(get_file_pointer(), elem_map.data());
      if (error < 0) {
        // Clear out the vector...
        Ioss::MapContainer().swap(elemMap.map());
        pamgen_error(get_file_pointer(), __LINE__, myProcessor);
      }
      elemMap.set_map(elem_map.data(), elem_map.size(), 0, true);
    }
    else {
      unsupported("output element map");
    }
  }
  return elemMap;
}

void DatabaseIO::compute_block_membership__(Ioss::SideBlock *         sideblock,
                                            std::vector<std::string> &block_membership) const
{
  Ioss::IntVector block_ids(elementBlockCount);
  if (elementBlockCount == 1) {
    block_ids[0] = 1;
  }
  else {
    Ioss::IntVector element_side;
    sideblock->get_field_data("element_side", element_side);

    size_t              number_sides = element_side.size() / 2;
    Ioss::ElementBlock *block        = nullptr;
    for (size_t iel = 0; iel < number_sides; iel++) {
      int elem_id = element_side[2 * iel]; // Vector contains both element and side.
      elem_id     = element_global_to_local(elem_id);
      if (block == nullptr || !block->contains(elem_id)) {
        block = get_region()->get_element_block(elem_id);
        assert(block != nullptr);
        int block_order        = block->get_property("original_block_order").get_int();
        block_ids[block_order] = 1;
      }
    }
  }

  // Synchronize among all processors....
  if (isParallel) {
    util().global_array_minmax(block_ids, Ioss::ParallelUtils::DO_MAX);
  }

  const Ioss::ElementBlockContainer &element_blocks = get_region()->get_element_blocks();

  for (int i = 0; i < elementBlockCount; i++) {
    if (block_ids[i] == 1) {
      Ioss::ElementBlock *block = element_blocks[i];
      block_membership.push_back(block->name());
    }
  }
}

int DatabaseIO::get_side_connectivity(const Ioss::SideBlock *fb, int id, int, int *fconnect,
                                      size_t /* data_size */) const
{
  // Get size of data stored on the file...
  int ierr;
  int number_sides;
  int number_distribution_factors;
  ierr =
      im_ex_get_side_set_param(get_file_pointer(), id, &number_sides, &number_distribution_factors);
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

  Ioss::IntVector is_valid_side;
  Ioss::Utils::calculate_sideblock_membership(is_valid_side, fb, 4, element.data(), side.data(),
                                              number_sides, get_region());

  Ioss::IntVector     elconnect;
  int                 elconsize  = 0;       // Size of currently allocated connectivity block
  Ioss::ElementBlock *conn_block = nullptr; // Block that we currently
  // have connectivity for

  Ioss::ElementBlock *block = nullptr;

  Ioss::IntVector side_elem_map; // Maps the side into the elements
  // connectivity array
  int current_side = -1;
  int nelnode      = 0;
  int nfnodes      = 0;
  int ieb          = 0;
  int offset       = 0;
  for (int iel = 0; iel < number_sides; iel++) {
    if (is_valid_side[iel] == 1) {

      int elem_id = element[iel];

      // ensure we have correct connectivity
      block = get_region()->get_element_block(elem_id);
      assert(block != nullptr);
      assert(block->topology() != nullptr);

      if (conn_block != block) {
        int nelem = block->entity_count();
        nelnode   = block->topology()->number_nodes();
        // Used to map element number into position in connectivity array.
        // E.g., element 97 is the (97-offset)th element in this block and
        // is stored in array index (97-offset-1).
        offset = block->get_offset() + 1;
        if (elconsize < nelem * nelnode) {
          elconsize = nelem * nelnode;
          elconnect.resize(elconsize);
        }
        get_field_internal(block, block->get_field("connectivity"), &elconnect[0],
                           nelem * nelnode * sizeof(int));
        conn_block   = block;
        current_side = -1;
      }

      // NOTE: Element connectivity is returned with nodes in global id space
      if (current_side != side[iel]) {
        int side_id   = side[iel];
        side_elem_map = block->topology()->boundary_connectivity(side_id);
        current_side  = side[iel];
        assert(block->topology()->boundary_type(side[iel]) != nullptr);
        nfnodes = block->topology()->boundary_type(side[iel])->number_nodes();
      }
      for (int inode = 0; inode < nfnodes; inode++) {
        int global_node = elconnect[(elem_id - offset) * nelnode + side_elem_map[inode]];
        fconnect[ieb++] = global_node;
      }
    }
  }
  return ierr;
}

// Get distribution factors for the specified side block
int DatabaseIO::get_side_distributions(const Ioss::SideBlock *fb, int id, int my_side_count,
                                       double *dist_fact, size_t /* data_size */) const
{
  // NOTE: In pamgen, all distribution factors are 1.0
  const Ioss::ElementTopology *ftopo   = fb->topology();
  int                          nfnodes = ftopo->number_nodes();
  for (int i = 0; i < nfnodes * my_side_count; i++)
    dist_fact[i] = 1.0;
  return 0;
}

//------------------------------------------------------------------------
} // namespace Iopg
namespace {
  void separate_surface_element_sides(Ioss::IntVector &element, Ioss::IntVector &sides,
                                      Ioss::Region *region, Iopg::TopologyMap &topo_map,
                                      Iopg::TopologyMap &    side_map,
                                      Ioss::SurfaceSplitType split_type)
  {
    if (!element.empty()) {
      Ioss::ElementBlock *block = nullptr;
      // Topology of sides in current element block
      const Ioss::ElementTopology *common_ftopo = nullptr;
      const Ioss::ElementTopology *topo         = nullptr; // Topology of current side
      int                          current_side = -1;

      for (size_t iel = 0; iel < element.size(); iel++) {
        int elem_id = element[iel];
        if (block == nullptr || !block->contains(elem_id)) {
          block = region->get_element_block(elem_id);
          assert(block != nullptr);

          // nullptr if hetero sides on element
          common_ftopo = block->topology()->boundary_type(0);
          if (common_ftopo != nullptr)
            topo = common_ftopo;
          current_side = -1;
        }

        if (common_ftopo == nullptr && sides[iel] != current_side) {
          current_side = sides[iel];
          topo         = block->topology()->boundary_type(sides[iel]);
          assert(topo != nullptr);
        }
        std::pair<std::string, const Ioss::ElementTopology *> name_topo;
        if (split_type == Ioss::SPLIT_BY_TOPOLOGIES) {
          name_topo = std::make_pair(block->topology()->name(), topo);
        }
        else if (split_type == Ioss::SPLIT_BY_ELEMENT_BLOCK) {
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
} // namespace
