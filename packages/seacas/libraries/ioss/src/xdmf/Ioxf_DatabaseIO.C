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

#include <xdmf/Ioxf_DatabaseIO.h>

#include <string>
#include <cstring>
#include <algorithm>
#include <vector>
#include <map>
#include <iterator>

#include <string>

#include <Ioss_CodeTypes.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_SerializeIO.h>
#include <Ioss_ElementTopology.h>

#include <assert.h>

// The following eliminates the use of the "__ FILE __" variable which
// was causing code bloat on janus.
static const char *Version = "2009/08/18";

// ========================================================================
// Static internal helper functions
// ========================================================================
namespace {
  bool set_id(const Ioss::GroupingEntity *entity, const char *type, Ioxf::EntityIdSet *idset);
  int  get_id(const Ioss::GroupingEntity *entity, char type, Ioxf::EntityIdSet *idset);
  void xdmf_error(const std::string &msg, int lineno, int processor) {
    std::ostringstream errmsg;
    errmsg << "XDMF error " << "processor " << processor << " in Ioxf:: line " << lineno <<
	      "\n    file/msg: " << msg << " in file '" << Version << "'" << "\n";
    IOSS_ERROR(errmsg);
  }
  
  void output_only()
  {
    std::ostringstream errmsg;
    errmsg << "XDMF Error: XDMF can only be used for output databases; "
	   << "however, a function used for database *input* was called which is not "
	   << "allowed.\n";
    IOSS_ERROR(errmsg);
  }

  // Determines whether the input map is sequential (map[i] == i)
  // Assumes that map is '1-based', size stored in [0]
  bool find_displacement_field(Ioss::NameList &fields,
			       const Ioss::GroupingEntity *block,
			       int ndim, std::string *disp_name);

  int field_warning(const Ioss::GroupingEntity *ge,
		    const Ioss::Field &field, const std::string& inout);

  void clean_out(const Ioxf::XmlContainer &container)
  {
    for (size_t i=0; i < container.size(); i++)
      delete container[i];
  }

#ifndef NDEBUG
  bool check_block_order(const Ioss::ElementBlockContainer &blocks);
#endif
}

namespace Ioxf {
  // ========================================================================
  // Internal typedefs/structs/classes used in some algorithms.
  // ========================================================================
  typedef std::map<std::pair<const Ioss::ElementTopology*,
			     const Ioss::ElementTopology*>, int,
		   std::less<std::pair<const Ioss::ElementTopology*,
				       const Ioss::ElementTopology*> > > TopologyMap;
  typedef TopologyMap::value_type TopoMapPair;

  // ========================================================================
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
			 Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
			 const Ioss::PropertyManager &properties) :
    Ioss::DatabaseIO(region, filename, db_usage, communicator, properties),
    databaseTitle(""), spatialDimension(0),
    nodeCount(0), elementCount(0),
    nodeBlockCount(0), elementBlockCount(0), nodesetCount(0), sidesetCount(0),
    nodeCmapIds(NULL), nodeCmapNodeCnts(NULL),
    elemCmapIds(NULL), elemCmapElemCnts(NULL), commsetNodeCount(0),
    commsetElemCount(0),
    elementTruthTable(NULL), nodesetTruthTable(NULL), sidesetTruthTable(NULL),
    sequentialNG2L(true), sequentialEG2L(true), fileExists(false)
  {
    // A history file is only written on processor 0...
    if (db_usage == Ioss::WRITE_HISTORY)
      isParallel = false;

    dbState = Ioss::STATE_UNKNOWN;

    // Initalize XDMF arrays
    Hdf = new XdmfHDF();
    MainXML =  new std::ostringstream();
    NumOfIterations = 0;

    // Construct the HDF and XML filenames for XDMF
    std::string decoded_name = util().decode_filename(filename, isParallel);

    hdfname = Ioss::FileInfo(decoded_name+".h5");
    xmlname = Ioss::FileInfo(decoded_name+".xmf");
  }


  DatabaseIO::~DatabaseIO()
  {
    delete [] elementTruthTable;
    delete [] nodesetTruthTable;
    delete [] sidesetTruthTable;
    delete [] nodeCmapIds;
    delete [] nodeCmapNodeCnts;
    delete [] elemCmapIds;
    delete [] elemCmapElemCnts;

    clean_out(BlockGridXmls);
    clean_out(BlockParameterXmls);
    clean_out(BlockXmls);
    clean_out(BlockElementVarXmls);
    clean_out(BlockExtraAttributeXmls);
    clean_out(BlockNodeVarXmls);
    clean_out(BlockFinishXmls);

    delete Hdf;
    delete MainXML;
  }

  void DatabaseIO::finalize()
  {
    H5close();
  }


  void DatabaseIO::WriteHdfDataset(const std::string &FinalName, XdmfArray *ScalarArray, int lineno) const
  {

    Hdf->CopyType( ScalarArray );
    Hdf->CopyShape( ScalarArray );
    Hdf->SetUseSerialFile(1);
    Hdf->Open( FinalName.c_str(), "rw" );
    
    if ( Hdf->CreateDataset( FinalName.c_str() ) == XDMF_FAIL) {
      std::string errmsg = FinalName;
      errmsg = "could not Create DataSet " + errmsg;
      xdmf_error(errmsg, lineno, myProcessor);
    }
    if ( Hdf->Write( ScalarArray ) == NULL) {
      std::string errmsg = FinalName;
      errmsg = "could not Write DataSet " + errmsg;
      xdmf_error(errmsg, lineno, myProcessor);
    }
    Hdf->Close();
  }

  void DatabaseIO::InitXML(std::ostringstream *XML)
  {
    if( XML){
      // Title...
      Ioss::Region *region = get_region();
      std::string title_str;
      if (region->property_exists("title")) {
	title_str = region->get_property("title").get_string();
      } else {
	title_str = "Sierra Output Default Title";
      }
      
      *XML << "<?xml version=\"1.0\" ?>" << endl;
      *XML << "<Xdmf>" << endl;
      *XML << " <Domain>" << endl << endl;
      *XML << "<!-- Raw Information from SierraFrameWork -->" << endl;
      *XML << "<?ExodusII Title=\"" << title_str << "\"" << endl;
      *XML << "\tElements=\"" << elementCount << "\"" << endl;
      *XML << "\tElementBlocks=\"" << elementBlockCount << "\"" << endl;
      *XML << "\tNodes=\"" << nodeCount << "\"" << endl;
      *XML << "\tNodeSets=\"" << nodesetCount << "\"" << endl;
      *XML << "\tSideSets=\"" << sidesetCount << "\"" << endl;
      *XML << "\tspatialDimension=\"" << spatialDimension << "\"" << endl;
      if (nodeCount > 0) {
	*XML << "\tcoordXName=\"x\"" << endl;
	*XML << "\tcoordYName=\"y\"" << endl;
	if (spatialDimension == 3)
	  *XML << "\tcoordZName=\"z\"" << endl;
      }
      *XML << "\t?>" << endl << endl;
    }
  }
  
  void DatabaseIO::put_qa()
  {
  }

  void DatabaseIO::put_info()
  {
  }

  void DatabaseIO::read_meta_data()
  {
    output_only();
  }

  void DatabaseIO::read_region()
  {
    output_only();
  }

  void DatabaseIO::get_step_times()
  {
    output_only();
  }

  void DatabaseIO::read_communication_metadata()
  {
    output_only();
  }

  const Ioss::MapContainer& DatabaseIO::get_node_map() const
  {
    // Allocate space for node number map and read it in...
    // Can be called multiple times, allocate 1 time only

    if (nodeMap.empty()) {
      // Cast away the 'const'.  Conceptually, this is a const operation
      // since it doesn't change the clients view of the database.
      DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
      new_this->nodeMap.resize(nodeCount+1);

      if (is_input()) {
	output_only();
      } else {
	// Output database; nodeMap not set yet... Build a default map.
	for (int i=1; i < nodeCount+1; i++) {
	  new_this->nodeMap[i] = i;
	}
	// Sequential map
	new_this->nodeMap[0] = -1;
      }
    }
    return nodeMap;
  }

  const Ioss::MapContainer& DatabaseIO::get_element_map() const
  {
    // Allocate space for elemente number map and read it in...
    // Can be called multiple times, allocate 1 time only
    if (elementMap.empty()) {
      // Cast away the 'const'.  Conceptually, this is a const operation
      // since it doesn't change the clients view of the database.
      DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
      new_this->elementMap.resize(elementCount+1);

      if (is_input()) {
	output_only();
      } else {
	// Output database; elementMap not set yet... Build a default map.
	for (int i=1; i < elementCount+1; i++) {
	  new_this->elementMap[i] = i;
	}
	// Sequential map
	new_this->sequentialEG2L = true;
	new_this->elementMap[0] = -1;
      }
    }
    return elementMap;
  }

  void DatabaseIO::get_nodeblocks() const
  {
    output_only();
  }

  void DatabaseIO::get_elemblocks() const
  {
    output_only();
  }

  void DatabaseIO::get_sidesets() const
  {
    output_only();
  }

  void DatabaseIO::get_nodesets() const
  {
    output_only();
  }

  void DatabaseIO::get_commsets() const
  {
    output_only();
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::Region* /* region */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock* /* nb */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeBlock* /* nb */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::FaceBlock* /* nb */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock* /* eb */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet* /* ns */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeSet* /* ns */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::FaceSet* /* ns */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::ElementSet* /* ns */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::SideSet* /* fs */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::CommSet* /* cs */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock* /* fb */,
				     const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    output_only();
    return 0;
  }
  
  int DatabaseIO::get_side_connectivity(const Ioss::EntityBlock* /* fb */,
					int /* id */, int /* side_count */,
					int */* fconnect */,
					size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  // Get distribution factors for the specified side block
  int DatabaseIO::get_side_distributions(const Ioss::EntityBlock* /* fb */,
					 int /* id */, int /* side_count */,
					 double */* dist_fact */,
					 size_t /* data_size */) const
  {
    output_only();
    return 0;
  }

  //------------------------------------------------------------------------
  int64_t DatabaseIO::put_field_internal(const Ioss::Region* /* region */,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    // For now, assume that all TRANSIENT fields on a region
    // are REDUCTION fields (1 value).  We need to gather these
    // and output them all at one time.  The storage location is a
    // 'globalVariables' array
    {
      Ioss::SerializeIO	serializeIO__(this);

      Ioss::Field::RoleType role = field.get_role();
      int num_to_get = field.verify(data_size);

      if ((role == Ioss::Field::TRANSIENT || role == Ioss::Field::REDUCTION) &&
	  num_to_get == 1) {
	store_reduction_field("R", field, get_region(), data);
      }
      else {
	std::ostringstream errmsg;
	errmsg << "Can not handle non-TRANSIENT or non-REDUCTION fields on regions ";
	IOSS_ERROR(errmsg);
	//    num_to_get = field_warning(region, field, "output");
      }
      return num_to_get;
    }
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock* nb,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    {
      Ioss::SerializeIO	serializeIO__(this);

      int num_to_get = field.verify(data_size);
      if (num_to_get > 0) {

	Ioss::Field::RoleType role = field.get_role();

	if (role == Ioss::Field::MESH) {
	  if (field.get_name() == "mesh_model_coordinates") {
	    // Data required by upper classes store x0, y0, z0, ... xn, yn, zn
	    // Data stored is in XDMF the same format

	    // Cast 'data' to correct size -- double
	    double *rdata = (double*)data;

	    XdmfInt64 Dimensions[2];
	    Dimensions[0] = num_to_get;
	    Dimensions[1] = spatialDimension;

	    XdmfArray *ScalarArray = new XdmfArray();

	    if( sizeof(double) == 8 ){
	      ScalarArray->SetNumberType( XDMF_FLOAT64_TYPE );
	    } else {
	      ScalarArray->SetNumberType( XDMF_FLOAT32_TYPE );
	    }

	    ScalarArray->SetNumberOfElements( num_to_get );
	    ScalarArray->SetDataPointer(&rdata[0]);
	    ScalarArray->SetShape( 2, Dimensions ); // Rank, Dimensions

	    std::string FinalName(hdfname.tailname());
	    FinalName += ":/Geometry";
	    WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	    delete ScalarArray;



	  } else if (field.get_name() == "ids") {
	    // The ids coming in are the global ids; their position is the
	    // local id -1 (That is, data[0] contains the global id of local
	    // node 1)

	    // Another 'const-cast' since we are modifying the database just
	    // for efficiency; which the client does not see...
	    DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
	    new_this->handle_node_ids((int*)data, num_to_get);
	  } else if (field.get_name() == "connectivity") {
	    // Do nothing, just handles an idiosyncracy of the GroupingEntity
	  } else {
	    return field_warning(nb, field, "mesh output");
	  }

	} else if (role == Ioss::Field::TRANSIENT) {
	  // Check if the specified field exists on this node block.
	  // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	  // exist on the database as scalars with the appropriate
	  // extensions.

	  // Transfer each component of the variable into 'data' and then
	  // output.  Need temporary storage area of size 'number of
	  // nodes in this block.
	  write_nodal_transient_field("N", field, nb, num_to_get, data);

	} else if (role == Ioss::Field::REDUCTION) {
	  store_reduction_field("N", field, nb, data);
	}
      }
      return num_to_get;
    }
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock* fs,
				     const Ioss::Field& field,
				     void */* data */, size_t /* data_size */) const
  {
    int num_to_get = field_warning(fs, field, "output");
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock* fs,
				     const Ioss::Field& field,
				     void */* data */, size_t /* data_size */) const
  {
    int num_to_get = field_warning(fs, field, "output");
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock* eb,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    {
      Ioss::SerializeIO	serializeIO__(this);

      int num_to_get = field.verify(data_size);
      std::ostringstream *XML=NULL;
      std::string block_name = eb->name();

      if (num_to_get > 0) {

	// Get the element block id and element count
	get_id(eb, 'E', &ids_);
	int element_count = eb->get_property("entity_count").get_int();
	Ioss::Field::RoleType role = field.get_role();

	if (role == Ioss::Field::MESH) {
	  // Handle the MESH fields required for an ExodusII file model.
	  // (The 'genesis' portion)
	  if (field.get_name() == "connectivity") {
	    if (element_count > 0) {
	      // Map element connectivity from global node id to local node id.
	      // Do it in 'data' ...
	      int* connect = (int*)data;
	      int element_nodes =
		eb->get_property("topology_node_count").get_int();

	      if (!sequentialNG2L) {
		assert(field.transformed_storage()->component_count() == element_nodes);

		for (int i=0; i < num_to_get * element_nodes; i++) {
		  int global_id = connect[i];
		  connect[i] = node_global_to_local(global_id, true);
		}
	      }

	      XdmfArray *ScalarArray = new XdmfArray();

	      ScalarArray->SetNumberType( XDMF_INT32_TYPE );
	      // ScalarArray->SetNumberOfElements( num_to_get );
	      ScalarArray->SetDataPointer(connect);

	      XdmfInt64 Dimensions[2];
	      Dimensions[0] = num_to_get;
	      Dimensions[1] = element_nodes;
	      ScalarArray->SetShape( 2, Dimensions );

	      std::string FinalName(hdfname.tailname());
	      FinalName += ":/" + block_name + "/Connections";
	      WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	      delete ScalarArray;

	    }
	  }

	  else if (field.get_name() == "thickness" || field.get_name() == "nodal_thickness" ||
		   field.get_name() == "radius" ||
		   field.get_name() == "area" || field.get_name() == "attribute") {
	    int attribute_count = eb->get_property("attribute_count").get_int();
	    if (element_count > 0 && attribute_count > 0) {

	      std::string name = field.get_name();

	      int ElementBlockIndex =  get_xml_stream(block_name);

	      XML = BlockGridXmls[ElementBlockIndex];

	      *XML << "\t<Attribute Name=\"" << name << "\"" << endl;
	      *XML << "\t\tCenter=\"" << "Grid" << "\"" << endl;
	      *XML << "\t\tType=\"" << "Scalar" << "\">" << endl;
	      *XML << "\t\t<DataStructure" << endl;
	      *XML << "\t\t\tFormat=\"HDF\"" << endl;
	      *XML << "\t\t\tDimensions=\"" << num_to_get << "\"" << endl;
	      *XML << "\t\t\tDataType=\""<< "Float" << "\" >" << endl;
	      *XML << "\t\t" << hdfname.tailname()
		   << ":/" << block_name << "/" << name << endl;
	      *XML << "\t\t</DataStructure>" << endl;
	      *XML << "\t</Attribute>" << endl;


	      XdmfArray *ScalarArray = new XdmfArray();

	      if( sizeof(double) == 8 ){
		ScalarArray->SetNumberType( XDMF_FLOAT64_TYPE );
	      } else {
		ScalarArray->SetNumberType( XDMF_FLOAT32_TYPE );
	      }

	      ScalarArray->SetNumberOfElements( num_to_get );
	      ScalarArray->SetDataPointer(data);

	      std::string FinalName(hdfname.tailname());
	      FinalName += ":/" + block_name + "/" + name;
	      WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	      delete ScalarArray;
	    }

	  } else if (field.get_name() == "ids") {
	    // Another 'const-cast' since we are modifying the database just
	    // for efficiency; which the client does not see...
	    DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
	    new_this->handle_element_ids(eb, (int*)data, num_to_get);

	  } else {
	    IOSS_WARNING << " ElementBlock "
			 << eb->name()
			 << ". Unknown field " << field.get_name();
	    num_to_get = -4; //IOSS_UNKNOWN_FIELD;
	  }
	} else if (role == Ioss::Field::TRANSIENT) {
	  // Check if the specified field exists on this element block.
	  // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	  // exist on the database as scalars with the appropriate
	  // extensions.

	  // Transfer each component of the variable into 'data' and then
	  // output.  Need temporary storage area of size 'number of
	  // elements in this block.
	  write_entity_transient_field("E", field, eb, element_count, data);

	} else if (role == Ioss::Field::REDUCTION) {
	  store_reduction_field("E", field, eb, data);
	}
      }
      return num_to_get;
    }
  }

  int DatabaseIO::handle_node_ids(int* ids, int num_to_get)
  {
    /*!
     * There are two modes we need to support in this routine:
     * 1. Initial definition of node map (local->global) and
     * reverseNodeMap (global->local).
     * 2. Redefinition of node map via 'reordering' of the original
     * map when the nodes on this processor are the same, but their
     * order is changed (or count because of ghosting)
     *
     * So, there will be two maps the 'nodeMap' map is a 'direct lookup'
     * map which maps current local position to global id and the
     * 'reverseNodeMap' is an associative lookup which maps the
     * global id to 'original local'.  There is also a
     * 'reorderNodeMap' which is direct lookup and maps current local
     * position to original local.

     * The ids coming in are the global ids; their position is the
     * "local id-1" (That is, data[0] contains the global id of local
     * node 1 in this node block).
     *
     * int local_position = reverseNodeMap[NodeMap[i+1]]
     * (the nodeMap and reverseNodeMap are 1-based)
     *
     * To determine which map to update on a call to this function, we
     * use the following hueristics:
     * -- If the database state is 'STATE_MODEL:', then update the
     *    'reverseNodeMap' and 'nodeMap'
     *
     * -- If the database state is not STATE_MODEL, then leave the
     *    'reverseNodeMap' and 'nodeMap' alone since they correspond to the
     *    information already written to the database. [May want to add a
     *    STATE_REDEFINE_MODEL]
     *
     * -- In both cases, update the reorderNodeMap
     *
     * NOTE: The mapping is done on TRANSIENT fields only; MODEL fields
     *       should be in the orginal order...
     */

    // The ids coming in are the global ids; their position is the
    // local id -1 (That is, data[0] contains the global id of local
    // node 1)

    if (dbState == Ioss::STATE_MODEL) {
      // This optimization only valid if entire field available
      if (sequentialNG2L && num_to_get == nodeCount) {
	for (int i=0; i < num_to_get; i++) {
	  if (i+1 != ids[i]) {
	    sequentialNG2L = false;
	    break;
	  }
	}
      }

      if (!sequentialNG2L || num_to_get != nodeCount) {
	sequentialNG2L = false;
	Ioss::Map::build_reverse_map(&reverseNodeMap, ids, num_to_get, 0,
				     "node", myProcessor);
      }

      // Only a single nodeblock and all set
      if (num_to_get == nodeCount) {
	assert(sequentialNG2L || (int)reverseNodeMap.size() == nodeCount);
      }
      assert(get_region()->get_property("node_block_count").get_int() == 1);

      XdmfArray *ScalarArray = new XdmfArray();

      ScalarArray->SetNumberType( XDMF_INT32_TYPE );
      ScalarArray->SetNumberOfElements( num_to_get );
      ScalarArray->SetDataPointer(ids);

      std::string FinalName(hdfname.tailname());
      FinalName += ":/NodeData/ids";

      WriteHdfDataset(FinalName, ScalarArray, __LINE__);
      delete ScalarArray;

    }

    build_node_reorder_map(ids, num_to_get);
    return num_to_get;
  }

  int DatabaseIO::handle_element_ids(const Ioss::ElementBlock *eb,
				     int* ids, int num_to_get)
  {
    /*!
     * There are two modes we need to support in this routine:
     * 1. Initial definition of element map (local->global) and
     * reverseElementMap (global->local).
     * 2. Redefinition of element map via 'reordering' of the original
     * map when the elements on this processor are the same, but their
     * order is changed.
     *
     * So, there will be two maps the 'elementMap' map is a 'direct lookup'
     * map which maps current local position to global id and the
     * 'reverseElementMap' is an associative lookup which maps the
     * global id to 'original local'.  There is also a
     * 'reorderElementMap' which is direct lookup and maps current local
     * position to original local.

     * The ids coming in are the global ids; their position is the
     * local id -1 (That is, data[0] contains the global id of local
     * element 1 in this element block).  The 'model-local' id is
     * given by eb_offset + 1 + position:
     *
     * int local_position = reverseElementMap[ElementMap[i+1]]
     * (the elementMap and reverseElementMap are 1-based)
     *
     * But, this assumes 1..numel elements are being output at the same
     * time; we are actually outputting a blocks worth of elements at a
     * time, so we need to consider the block offsets.
     * So... local-in-block position 'i' is index 'eb_offset+i' in
     * 'elementMap' and the 'local_position' within the element
     * blocks data arrays is 'local_position-eb_offset'.  With this, the
     * position within the data array of this element block is:
     *
     * int eb_position =
     * reverseElementMap[elementMap[eb_offset+i+1]]-eb_offset-1
     *
     * To determine which map to update on a call to this function, we
     * use the following hueristics:
     * -- If the database state is 'Ioss::STATE_MODEL:', then update the
     *    'reverseElementMap'.
     * -- If the database state is not Ioss::STATE_MODEL, then leave
     *    the 'reverseElementMap' alone since it corresponds to the
     *    information already written to the database. [May want to add
     *    a Ioss::STATE_REDEFINE_MODEL]
     * -- Always update elementMap to match the passed in 'ids'
     *    array.
     *
     * NOTE: the maps are built an element block at a time...
     * NOTE: The mapping is done on TRANSIENT fields only; MODEL fields
     *       should be in the orginal order...
     */

    // Overwrite this portion of the 'elementMap', but keep other
    // parts as they were.  We are adding elements starting at position
    // 'eb_offset+offset' and ending at
    // 'eb_offset+offset+num_to_get'. If the entire block is being
    // processed, this reduces to the range 'eb_offset..eb_offset+element_count'
    if (elementMap.empty()) {
      elementMap.resize(elementCount+1);
      sequentialEG2L = true;
    }

    assert(static_cast<int>(elementMap.size()) == elementCount+1);
    elementMap[0] = -1;

    int eb_offset = eb->get_offset();

    for (int i=0; i < num_to_get; i++) {
      int local_id = eb_offset + i + 1;
      elementMap[local_id] = ids[i];
      if (local_id != ids[i]) {
	sequentialEG2L = false;
	elementMap[0] = 1;
      }
    }

    // Now, if the state is Ioss::STATE_MODEL, update the reverseElementMap
    if (dbState == Ioss::STATE_MODEL) {
      Ioss::Map::build_reverse_map(&reverseElementMap, ids, num_to_get,
				   eb_offset, "element", myProcessor);

      // Strongest assertion we can make is that size of map <=
      // elementCount
      assert((int)reverseElementMap.size() <= elementCount);

      // Output this portion of the element number map
      std::ostringstream *XML=NULL;
      std::string block_name = eb->name();

      int ElementBlockIndex =  get_xml_stream(block_name);

      XML = BlockNodeVarXmls[ElementBlockIndex];
      *XML << "\t<Attribute Name=\"" << "ids" << "\"" << endl;
      *XML << "\t\tCenter=\"" << "Cell" << "\"" << endl;
      *XML << "\t\tType=\"" << "Scalar" << "\">" << endl;               *XML << "\t\t<DataStructure" << endl;
      *XML << "\t\t\tFormat=\"HDF\"" << endl;
      *XML << "\t\t\tDimensions=\"" << num_to_get << "\"" << endl;
      *XML << "\t\t\tDataType=\""<< "Int" << "\" >" << endl;
      *XML << "\t\t" << hdfname.tailname()
	   << ":/" << block_name << "/" << "ids" << endl;
      *XML << "\t\t</DataStructure>" << endl;
      *XML << "\t</Attribute>" << endl;


      XdmfArray *ScalarArray = new XdmfArray();

      ScalarArray->SetNumberType( XDMF_INT32_TYPE );
      ScalarArray->SetNumberOfElements( num_to_get );
      ScalarArray->SetDataPointer(ids);

      std::string FinalName(hdfname.tailname());
      FinalName += ":/" + block_name + "/ids";
      WriteHdfDataset(FinalName, ScalarArray, __LINE__);
      delete ScalarArray;

    }
    // Build the reorderElementMap which does a direct mapping from
    // the current topologies local order to the local order stored in
    // the database...  This is 0-based and used for remapping output
    // TRANSIENT fields. (Will also need one on input once we read fields)
    build_element_reorder_map(eb_offset, num_to_get);
    return num_to_get;
  }

  void DatabaseIO::write_nodal_transient_field(const char */* type */,
					       const Ioss::Field &field,
					       const Ioss::NodeBlock */* ge */,
					       int count,
					       void *data) const
  {
    Ioss::Field::BasicType ioss_type = field.get_type();
    assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER);
    double *rvar = (double*)data;
    int  *ivar = (int*)data;

    const Ioss::VariableType *var_type = field.transformed_storage();
    std::vector<double> temp(count);

    int step = get_region()->get_property("current_state").get_int();
    step = get_database_step(step);

    // get number of components, cycle through each component
    // and add suffix to base 'field_name'.  Look up index
    // of this name in 'nodalVariables' map
    int comp_count = var_type->component_count();
    for (int i=0; i < comp_count; i++) {
      std::string var_name = var_type->label_name(field.get_name(), i+1, get_field_separator());
      DatabaseIO *new_this = const_cast<DatabaseIO*>(this);

      new_this->nodalVariables[var_name];

      // Transfer from 'variables' array.
      int k = 0;
      int num_out = 0;
      for (int j=i; j < count*comp_count; j+=comp_count) {
	int where = reorderNodeMap[k++];
	if (where >= 0) {
	  assert(where < count);
	  if (ioss_type == Ioss::Field::REAL)
	    temp[where] = rvar[j];
	  else
	    temp[where] = ivar[j];
	  num_out++;
	}
      }
      assert(k == count);
      assert(num_out == nodeCount);

      // Write the variable...

      // Write a NodeBlock for index 0 for every time step
      std::ostringstream *XML=NULL;

      if (step == 1) {
	int num_elemets_blocks = BlockNodeVarXmls.size();
	for (int ii = 0 ; ii < num_elemets_blocks ; ii++ ){
	  XML = BlockNodeVarXmls[ii];
	  *XML << "\t<Attribute Name=\"" << var_name << "\"" << endl;
	  *XML << "\t\tCenter=\"" << "Node" << "\"" << endl;
	  *XML << "\t\tType=\"" << "Scalar" << "\">" << endl;
	  *XML << "\t\t<DataStructure" << endl;
	  *XML << "\t\t\tFormat=\"HDF\"" << endl;
	  *XML << "\t\t\tDimensions=\"" << num_out << "\"" << endl;
	  *XML << "\t\t\tDataType=\""<< "Float" << "\">" << endl;
	  *XML << "\t\t" << hdfname.tailname()
	       << ":/NodeData/Iteration/" << var_name << endl;

	  *XML << "\t\t</DataStructure>" << endl;
	  *XML << "\t</Attribute>" << endl;
	}
      }

      XdmfArray *ScalarArray = new XdmfArray();

      if( sizeof(double) == 8 ){
	ScalarArray->SetNumberType( XDMF_FLOAT64_TYPE );
      } else {
	ScalarArray->SetNumberType( XDMF_FLOAT32_TYPE );
      }

      ScalarArray->SetNumberOfElements(num_out);
      ScalarArray->SetDataPointer(&temp[0]);

      std::string FinalName(hdfname.tailname());
      FinalName += ":/NodeData/Iteration ";
      FinalName += Ioss::Utils::to_string(step);
      FinalName += "/";
      FinalName += var_name;
      WriteHdfDataset(FinalName, ScalarArray, __LINE__);
      delete ScalarArray;

    }

    if (step == 1 && (field.get_name().c_str()[0] == 'd' ||
		      field.get_name().c_str()[0] == 'D') &&
	comp_count == spatialDimension) {
      int num_elemets_blocks = BlockNodeVarXmls.size();
      for (int i = 0 ; i < num_elemets_blocks ; i++ ){
	std::ostringstream *XML= BlockNodeVarXmls[i];
	*XML << "\t<Attribute Name=\"DISPLACEMENT\" Type=\"Vector\" Center=\"Node\">" << endl;
	*XML << "\t\t<DataTransform Type=\"Function\" Format=\"XML\" Function=\"join($0 , $1 ,  $2)\""
	     << " Dimensions=\"" << nodeCount << " 3\">" << endl;

	for (int j=0; j < comp_count; j++) {
	  std::string var_name = var_type->label_name(field.get_name(), j+1, get_field_separator());

	  *XML << "\t\t\t<DataStructure Format=\"HDF\" Dimensions=\"" << nodeCount << "\" DataType=\"Float\">" << endl;
	  *XML << "\t\t\t\t" << hdfname.tailname()
	       << ":/NodeData/Iteration/" << var_name << endl;
	  *XML << "\t\t\t</DataStructure>" << endl;
	}

	*XML << "\t\t</DataTransform>" << endl;
	*XML << "\t</Attribute>" << endl;
      }
    }
  }

  void DatabaseIO::write_entity_transient_field(const char *type,
						const Ioss::Field& field,
						const Ioss::GroupingEntity *ge,
						int count,
						void *data) const
  {
    assert(type[0] == 'E' || type[0] == 'M' || type[0] == 'S');

    const Ioss::VariableType *var_type = field.transformed_storage();

    std::vector<double> temp(count);

    int step = get_region()->get_property("current_state").get_int();
    step = get_database_step(step);

    int eb_offset = 0;
    if (ge->type() == Ioss::ELEMENTBLOCK) {
      const Ioss::ElementBlock *elb = dynamic_cast<const Ioss::ElementBlock*>(ge);
      assert(elb != NULL);
      eb_offset = elb->get_offset();
    }

    Ioss::Field::BasicType ioss_type = field.get_type();
    assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER);
    double *rvar = (double*)data;
    int  *ivar = (int*)data;

    // get number of components, cycle through each component
    // and add suffix to base 'field_name'.  Look up index
    // of this name in 'elementVariables' map
    int comp_count = var_type->component_count();
    for (int i=0; i < comp_count; i++) {
      std::string var_name = var_type->label_name(field.get_name(), i+1, get_field_separator());

      // Transfer from 'variables' array.  Note that the
      // 'reorderElementMap has '1..numel' ids in it, but the 'temp'
      // array is stored in 'element block local id space', so we need
      // to add/subtract the element block offset to keep things
      // straight.
      int k = eb_offset;
      int where = 0;
      for (int j=i; j < count*comp_count; j+=comp_count) {
	// Map to storage location.

	if (type[0] == 'E')
	  where = reorderElementMap[k++] - eb_offset;
	else
	  where = k++;

	assert(where >= 0 && where < count);
	if (ioss_type == Ioss::Field::REAL)
	  temp[where] = rvar[j];
	else
	  temp[where] = ivar[j];
      }
      assert(k-eb_offset == count);

      // Write the variable...
      int id;
      id  = get_id(ge, type[0], &ids_);

      if (type[0] == 'E') {
	std::string block_name = ge->name();

	std::string FinalName(hdfname.tailname().c_str());
	FinalName += ":/" + block_name + "/Iteration " + Ioss::Utils::to_string(step) + "/" + var_name;
	XdmfArray *ScalarArray = new XdmfArray();
	if( sizeof(double) == 8 )
	  ScalarArray->SetNumberType( XDMF_FLOAT64_TYPE );
	else
	  ScalarArray->SetNumberType( XDMF_FLOAT32_TYPE );

	ScalarArray->SetNumberOfElements(count);
	ScalarArray->SetDataPointer(&temp[0]);

	WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	delete ScalarArray;
      }
      //  Erik-- Still need to implement  'M'
      else if (type[0] == 'M') {
	XdmfArray *ScalarArray = new XdmfArray();
	if( sizeof(double) == 8 )
	  ScalarArray->SetNumberType( XDMF_FLOAT64_TYPE );
	else
	  ScalarArray->SetNumberType( XDMF_FLOAT32_TYPE );
	ScalarArray->SetNumberOfElements(count);
	ScalarArray->SetDataPointer(&temp[0]);

	std::string FinalName(hdfname.tailname().c_str());
	FinalName += ":/NodeSets/" + Ioss::Utils::to_string(id) + "/Iteration " + Ioss::Utils::to_string(step) + "/" + var_name;

	WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	delete ScalarArray;
      }
      else if (type[0] == 'S') {
	XdmfArray *ScalarArray = new XdmfArray();
	if( sizeof(double) == 8 )
	  ScalarArray->SetNumberType( XDMF_FLOAT64_TYPE );
	else
	  ScalarArray->SetNumberType( XDMF_FLOAT32_TYPE );
	ScalarArray->SetNumberOfElements(count);
	ScalarArray->SetDataPointer(&temp[0]);

	std::string FinalName(hdfname.tailname().c_str());
	FinalName += ":/SideSets/" + Ioss::Utils::to_string(id) + "/Iteration " + Ioss::Utils::to_string(step) + "/" + var_name;
	WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	delete ScalarArray;
      }
    }
  }

  void DatabaseIO::store_reduction_field(const char *type,
					 const Ioss::Field& field,
					 const Ioss::GroupingEntity *ge,
					 void *data) const
  {
    const Ioss::VariableType *var_type = field.transformed_storage();

    Ioss::Field::BasicType ioss_type = field.get_type();
    assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER);
    double *rvar = (double*)data;
    int  *ivar = (int*)data;

    // get number of components, cycle through each component
    // and add suffix to base 'field_name'.  Look up index
    // of this name in 'globalVariables' map
    int comp_count = var_type->component_count();
    DatabaseIO *new_this = const_cast<DatabaseIO*>(this);

    for (int i=0; i < comp_count; i++) {
      std::string var_name = var_type->label_name(field.get_name(), i+1, get_field_separator());

      // If this is an element block, need to prepend the block name
      // to avoid name collisions... May also need this for nodeblocks
      if (type[0] == 'E') {
	std::string ge_name = ge->name();
	var_name = ge_name + ":" + var_name;
      }
      int var_index = new_this->globalVariables[var_name];

      assert(var_index > 0);
      if ((int)globalValues.size() < var_index)
	globalValues.resize(new_this->globalVariables.size());
      assert((int)globalValues.size() >= var_index);


      // Transfer from 'variables' array.
      if (ioss_type == Ioss::Field::REAL)
	new_this->globalValues[var_index-1] = rvar[i];
      else
	new_this->globalValues[var_index-1] = ivar[i];
    }
  }

  void DatabaseIO::get_reduction_field(const char *type,
				       const Ioss::Field& field,
				       const Ioss::GroupingEntity */* ge */,
				       int /* count */, void *variables) const
  {
    assert(type[0] == 'R'); // Only region at this time.
    const Ioss::VariableType *var_type = field.raw_storage();

    Ioss::Field::BasicType ioss_type = field.get_type();
    assert(ioss_type == Ioss::Field::REAL || ioss_type == Ioss::Field::INTEGER);
    double *rvar = static_cast<double*>(variables);
    int  *ivar = static_cast<int*>(variables);

    // get number of components, cycle through each component
    // and add suffix to base 'field_name'.  Look up index
    // of this name in 'globalVariables' map
    int comp_count = var_type->component_count();
    for (int i=0; i < comp_count; i++) {
      std::string var_name = var_type->label_name(field.get_name(), i+1, get_field_separator());

      DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
      int var_index = new_this->globalVariables[var_name];

      assert(var_index > 0);
      assert((int)globalValues.size() >= var_index);

      // Transfer to 'variables' array.
      if (ioss_type == Ioss::Field::REAL)
	rvar[i] = globalValues[var_index-1];
      else
	ivar[i] = static_cast<int>(globalValues[var_index-1]);
    }
  }

  void DatabaseIO::write_reduction_fields() const
  {
    // Not supported...
  }

  void DatabaseIO::read_reduction_fields() const
  {
    output_only();
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet* ns,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    {
      Ioss::SerializeIO	serializeIO__(this);

      int entity_count = ns->get_property("entity_count").get_int();
      int num_to_get = field.verify(data_size);
      if (num_to_get > 0) {

	int id = get_id(ns, 'N', &ids_);
	Ioss::Field::RoleType role = field.get_role();

	if (role == Ioss::Field::MESH) {

	  if (field.get_name() == "ids") {
	    // Map node id from global node id to local node id.
	    // Do it in 'data' ...
	    int* ids = (int*)data;

	    if (!sequentialNG2L) {
	      for (int i=0; i < num_to_get; i++) {
		int global_id = ids[i];
		ids[i] = node_global_to_local(global_id, true);
	      }
	    }

	    XdmfArray *ScalarArray = new XdmfArray();
	    ScalarArray->SetNumberType( XDMF_INT32_TYPE );
	    ScalarArray->SetNumberOfElements( num_to_get );
	    ScalarArray->SetDataPointer(data);

	    std::string FinalName(hdfname.tailname());
	    FinalName += ":/NodeSets/" + Ioss::Utils::to_string(id) + "/Node";
	    WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	    delete ScalarArray;

	  } else if (field.get_name() == "distribution_factors") {

	    XdmfArray *ScalarArray = new XdmfArray();
	    if( sizeof(double) == 8 )
	      ScalarArray->SetNumberType( XDMF_FLOAT64_TYPE );
	    else
	      ScalarArray->SetNumberType( XDMF_FLOAT32_TYPE );
	    ScalarArray->SetNumberOfElements( num_to_get );
	    ScalarArray->SetDataPointer(data);

	    std::string FinalName(hdfname.tailname());
	    FinalName += ":/NodeSets/" + Ioss::Utils::to_string(id) + "/DistributionFactor";
	    WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	    delete ScalarArray;

	  } else {
	    num_to_get = field_warning(ns, field, "output");
	  }
	} else if (role == Ioss::Field::TRANSIENT) {
	  // Check if the specified field exists on this element block.
	  // Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	  // exist on the database as scalars with the appropriate
	  // extensions.

	  // Transfer each component of the variable into 'data' and then
	  // output.  Need temporary storage area of size 'number of
	  // elements in this block.
	  write_entity_transient_field("M", field, ns, entity_count, data);

	} else if (role == Ioss::Field::REDUCTION) {
	  store_reduction_field("M", field, ns, data);
	}
      }
      return num_to_get;
    }
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeSet* fs,
				     const Ioss::Field& field,
				     void */* data */, size_t /* data_size */) const
  {
    int num_to_get = field_warning(fs, field, "output");
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::FaceSet* fs,
				     const Ioss::Field& field,
				     void */* data */, size_t /* data_size */) const
  {
    int num_to_get = field_warning(fs, field, "output");
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementSet* fs,
				     const Ioss::Field& field,
				     void */* data */, size_t /* data_size */) const
  {
    int num_to_get = field_warning(fs, field, "output");
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::SideSet* fs,
				     const Ioss::Field& field,
				     void */* data */, size_t /* data_size */) const
  {
    int num_to_get = field_warning(fs, field, "output");
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::CommSet* cs,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    int num_to_get = field.verify(data_size);
    int entity_count = cs->get_property("entity_count").get_int();

    assert(num_to_get == entity_count);
    if (num_to_get == 0)
      return 0;

    // Return the <entity (node or side), processor> pair
    if (field.get_name() == "entity_processor") {

      // Check type -- node or side
      std::string type = cs->get_property("entity_type").get_string();

      // Allocate temporary storage space
      std::vector<int> entities(entity_count);
      std::vector<int> procs(entity_count);

      if (type == "node") {
	// Convert global node id to local node id and store in 'entities'
	int* entity_proc = (int*)data;
	int j=0;
	for (int i=0; i < entity_count; i++) {
	  int global_id = entity_proc[j++];
	  entities[i] = node_global_to_local(global_id, true);
	  procs[i] = entity_proc[j++];
	}
	/*
	  if (commsetNodeCount > 0) {
	  int ierr = ne_put_node_cmap(get_file_pointer(), get_id(cs, 'C', &ids_),
	  &entities[0], &procs[0], myProcessor);
	  if (ierr < 0)
	  xdmf_error(get_file_pointer(), __LINE__, myProcessor);
	  }
	*/

	if (commsetNodeCount == 1) {
	  // NOTE: The internal and border node maps must be output in one call.
	  //       In this routine, we only have one commset at a time and can't
	  //       construct the entire map at one time.  This is not really needed,
	  //       so for now we just skip if there is more than one commset.  If
	  //       this information is really needed, need to cache the information
	  //       until all commsets have been processed.  Also need to change
	  //       write_communication_metada() [Maybe, unless client sets correct
	  //       properties.]

	  // Construct the node map (internal vs. border).
	  // Border nodes are those in the communication map (use entities array)
	  // Internal nodes are the rest.  Allocate array to hold all nodes,
	  // initialize all to '1', then zero out the nodes in 'entities'.
	  // Iterate through array again and consolidate all '1's
	  std::vector<int> internal(nodeCount);
	  for (j=0; j < nodeCount; j++)
	    internal[j] = 1;
	  for (j=0; j < entity_count; j++)
	    internal[entities[j]-1] = 0;

	  int b = 0;
	  for (j=0; j < nodeCount; j++) {
	    if (internal[j] == 0) {
	      entities[b++] = j+1;
	    }
	  }

	  int k = 0;
	  for (j=0; j < nodeCount; j++) {
	    if (internal[j] == 1) {
	      internal[k++] = j+1;
	    }
	  }

#ifndef NDEBUG
	  int border_nodes   = 0;
	  int internal_nodes = 0;
	  if (get_region()->property_exists("border_node_count")) {
	    border_nodes = get_region()->get_property("border_node_count").get_int();
	    assert(border_nodes   == b);
	  }

	  if (get_region()->property_exists("internal_node_count")) {
	    internal_nodes = get_region()->
	      get_property("internal_node_count").get_int();
	    assert(internal_nodes == k);
	  }
#endif
	  /*
	    int ierr = ne_put_node_map(get_file_pointer(), &internal[0], &entities[0], NULL,
	    myProcessor);
	    if (ierr < 0)
	    xdmf_error(get_file_pointer(), __LINE__, myProcessor);
	  */
	}

      } else if (type == "side") {
	std::vector<int> sides(entity_count);
	int* entity_proc = (int*)data;
	int j=0;
	for (int i=0; i < entity_count; i++) {
	  // Assume klugy side id generation.
	  int global_id = entity_proc[j] / 10;
	  entities[i] = element_global_to_local(global_id);
	  sides[i] = entity_proc[j++] % 10;
	  procs[i] = entity_proc[j++];
	}

#if 0
	int offset = 0;
	for (int ics=0; ics < commsetElemCount; ics++) {
	  int ierr = ne_put_elem_cmap(get_file_pointer(), elemCmapIds[ics],
				      &entities[offset], &sides[offset],
				      &procs[offset], myProcessor);
	  if (ierr < 0)
	    xdmf_error(get_file_pointer(), __LINE__, myProcessor);
	  offset += elemCmapElemCnts[ics];
	}
#else
	/*
	  int ierr = ne_put_elem_cmap(get_file_pointer(), get_id(cs, 'C', &ids_),
	  &entities[0], &sides[0], &procs[0], myProcessor);
	  if (ierr < 0)
	  xdmf_error(get_file_pointer(), __LINE__, myProcessor);
	*/
#endif

	// Construct the element map (internal vs. border).
	// Border elements are those in the communication map (use entities array)
	// Internal elements are the rest.  Allocate array to hold all elements,
	// initialize all to '1', then zero out the elements in 'entities'.
	// Iterate through array again and consolidate all '1's
	std::vector<int> internal(elementCount);
	for (j=0; j < elementCount; j++)
	  internal[j] = 1;
	for (j=0; j < entity_count; j++)
	  internal[entities[j]-1] = 0;

	int b = 0;
	for (j=0; j < elementCount; j++) {
	  if (internal[j] == 0) {
	    entities[b++] = j+1;
	  }
	}

	int k = 0;
	for (j=0; j < elementCount; j++) {
	  if (internal[j] == 1) {
	    internal[k++] = j+1;
	  }
	}

#ifndef NDEBUG
	int border_elems   = 0;
	int internal_elems = 0;
	if (get_region()->property_exists("border_elem_count")) {
	  border_elems =
	    get_region()->get_property("border_elem_count").get_int();
	  assert(border_elems   == b);
	}

	if (get_region()->property_exists("internal_elem_count")) {
	  internal_elems =
	    get_region()->get_property("internal_elem_count").get_int();
	  assert(internal_elems == k);
	}
#endif
	/*
	  ierr = ne_put_elem_map(get_file_pointer(), &internal[0],
	  &entities[0], myProcessor);
	  if (ierr < 0)
	  xdmf_error(get_file_pointer(), __LINE__, myProcessor);
	*/

      }
    } else {
      num_to_get = field_warning(cs, field, "output");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock* fb,
				     const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    int num_to_get = field.verify(data_size);
    if (num_to_get > 0) {

      int id = get_id(fb, 'S', &ids_);

      int entity_count = fb->get_property("entity_count").get_int();
      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
	if (field.get_name() == "ids") {
	  // =============================================================
	  // NOTE: We have redundant ways of getting the data
	  // (element/side) out to the database.  The 'ids' field
	  // method relies on a numbering kluge, so for now trying the
	  // 'element_side' field...
	  // =============================================================
	} else if (field.get_name() == "distribution_factors") {
	  int node_count = fb->get_property("topology_node_count").get_int();
	  int df_to_get = node_count * num_to_get;

	  std::string FinalName(hdfname.tailname());
	  FinalName += ":/SideSets/" + Ioss::Utils::to_string(id) + "/DistributionFactor";

	  XdmfArray *ScalarArray = new XdmfArray();
	  if( sizeof(double) == 8 )
	    ScalarArray->SetNumberType( XDMF_FLOAT64_TYPE );
	  else
	    ScalarArray->SetNumberType( XDMF_FLOAT32_TYPE );
	  ScalarArray->SetNumberOfElements(df_to_get);
	  ScalarArray->SetDataPointer(data);

	  WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	  delete ScalarArray;

	} else if (field.get_name() == "element_side") {
#if 1
	  // In exodusII, the 'side block' is stored as a sideset.  A
	  // sideset has a list of elements and a corresponding local
	  // element side (1-based)

	  // The 'data' passed into the function is stored as a
	  // 2D vector e0,f0,e1,f1,... (e=element, f=side)

	  // To avoid overwriting the passed in data, we allocate
	  // two arrays to store the data for this sideset.

	  // The element_id passed in is the global id; we need to
	  // output the local id.

	  // Allocate space for local side number and element numbers
	  // numbers.
	  std::vector<int> element(num_to_get);
	  std::vector<int> side(num_to_get);
	  int *el_side = (int *)data;

	  int index = 0;
	  for (int i=0; i < num_to_get; i++) {
	    element[i] = element_global_to_local(el_side[index++]);
	    side[i]    = el_side[index++];
	  }


	  XdmfArray *ScalarArray = new XdmfArray();
	  ScalarArray->SetNumberType( XDMF_INT32_TYPE );
	  ScalarArray->SetNumberOfElements( num_to_get );
	  ScalarArray->SetDataPointer(&side[0]);


	  {
	    std::string FinalName(hdfname.tailname());
	    FinalName += ":/SideSets/" + Ioss::Utils::to_string(id) + "/Side";
	    WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	  }

	  delete ScalarArray;

	  {
	    ScalarArray = new XdmfArray();
	    ScalarArray->SetNumberType( XDMF_INT32_TYPE );
	    ScalarArray->SetNumberOfElements( num_to_get );
	    ScalarArray->SetDataPointer(&element[0]);

	    XdmfInt64        Dimensions[2];
	    Dimensions[0] = num_to_get;
	    Dimensions[1] = 1;
	    ScalarArray->SetShape( 2, Dimensions ); // Rank, Dimensions

	    std::string FinalName(hdfname.tailname());
	    FinalName += ":/SideSets/" + Ioss::Utils::to_string(id) + "/Element";
	    WriteHdfDataset(FinalName, ScalarArray, __LINE__);
	    delete ScalarArray;
	    // delete Hdf;
	  }

#endif
	} else if (field.get_name() == "connectivity") {
	  // Do nothing, just handles an idiosyncracy of the GroupingEntity
	} else {
	  num_to_get = field_warning(fb, field, "output");
	}
      } else if (role == Ioss::Field::TRANSIENT) {
	// Check if the specified field exists on this block.
	// Note that 'higher-order' storage types (e.g. SYM_TENSOR)
	// exist on the database as scalars with the appropriate
	// extensions.

	// Transfer each component of the variable into 'data' and then
	// output.  Need temporary storage area of size 'number of
	// entities in this block.
	write_entity_transient_field("S", field, fb, entity_count, data);

      } else if (role == Ioss::Field::REDUCTION) {
	store_reduction_field("S", field, fb, data);
      }
    }
    return num_to_get;
  }

  // ------------------------------------------------------------------------
  // Node and Element mapping functions.  The ExodusII database
  // stores ids in a local-id system (1..NUMNP), (1..NUMEL) but
  // Sierra wants entities in a global system. These routines
  // take care of the mapping from local <-> global

  int64_t DatabaseIO::node_local_to_global(int64_t local)  const
  {
    assert(local <= nodeCount && local > 0);
    const Ioss::MapContainer &node_map = get_node_map();
    int global = node_map[local];
    return global;
  }

  int64_t DatabaseIO::element_local_to_global(int64_t local)  const
  {
    assert(local <= elementCount && local > 0);
    const Ioss::MapContainer &element_map = get_element_map();
    int global = element_map[local];
    return global;
  }

  bool DatabaseIO::begin(Ioss::State state)
  {
    dbState = state;
    return true;
  }

  bool DatabaseIO::end(Ioss::State state)
  {
    // Transitioning out of state 'state'
    assert(state == dbState);
    switch (state) {
    case Ioss::STATE_DEFINE_MODEL:
      if (!is_input())
	write_meta_data();
      break;
    case Ioss::STATE_DEFINE_TRANSIENT:
      if (!is_input())
	write_results_metadata();
      break;
    default: // ignore everything else...
      if (!is_input()) {
	WriteXmlFile(NumOfIterations);
      }
      break;
    }

    {
      Ioss::SerializeIO	serializeIO__(this);

      dbState = Ioss::STATE_UNKNOWN;
    }

    return true;
  }

  // Default versions do nothing at this time...
  // Will be used for global variables...
  bool DatabaseIO::begin_state(Ioss::Region */* region */, int state, double /* time */)
  {
    Ioss::SerializeIO	serializeIO__(this);

    state = get_database_step(state);
    if (!is_input()) {
      // Zero global variable array...
      std::fill(globalValues.begin(), globalValues.end(), 0.0);

    } else {
      // Store reduction variables
      read_reduction_fields();
    }
    return true;
  }

  bool DatabaseIO::end_state(Ioss::Region */* region */, int state, double /* time */)
  {
    Ioss::SerializeIO	serializeIO__(this);

    if (!is_input()) {
      // write_reduction_fields();
    }
    if (!is_input()) {
      NumOfIterations = state;
      WriteXmlFile(NumOfIterations);
    }

    return true;
  }

  void DatabaseIO::write_meta_data()
  {
    Ioss::Region *region = get_region();

    // Get count of nodes, element blocks, nodesets, sidesets, ...
    std::vector<Ioxf::Block>   blocks;
    std::vector<Ioxf::NodeSet> nsets;
    std::vector<Ioxf::SideSet> ssets;

    {
      Ioss::SerializeIO	serializeIO__(this);

      {
	Ioss::NodeBlockContainer node_blocks = region->get_node_blocks();
	assert(node_blocks.size() == 1);
	spatialDimension = node_blocks[0]->
	  get_property("component_degree").get_int();
	nodeCount =        node_blocks[0]->
	  get_property("entity_count").get_int();
      }


      // Element Blocks --
      {
	Ioss::ElementBlockContainer element_blocks =
	  region->get_element_blocks();
	assert(check_block_order(element_blocks));
	Ioss::ElementBlockContainer::const_iterator I;
	// Set ids of all entities that have "id" property...
	for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
	  set_id(*I, "Element Block", &ids_);
	}

	elementCount = 0;
	for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
	  elementCount += (*I)->get_property("entity_count").get_int();
	  // Set ids of all entities that do not have "id" property...
	  get_id(*I, 'E', &ids_);
	  Ioxf::Block T(*(*I));
	  if (std::find(blocks.begin(), blocks.end(), T) == blocks.end()) {
	    blocks.push_back(T);
	  }
	}
	elementBlockCount = blocks.size();
      }
    }

    // SideSets ...
    Ioss::SideSetContainer sidesets = region->get_sidesets();
    Ioss::SideSetContainer::const_iterator I;
    // Set ids of all entities that have an existing "id" property


    // Get entity counts for all side blocks...
    std::vector<int> sset_entity_count;
    for (I=sidesets.begin(); I != sidesets.end(); ++I) {
      Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
      Ioss::SideBlockContainer::const_iterator J;

      for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
	int count = (*J)->get_property("entity_count").get_int();
	sset_entity_count.push_back(count);
      }
    }

    // Resolve count among all side blocks on all processors...
    // NOTE: This is a collective call.
    if (isParallel && !sset_entity_count.empty()) {
      util().global_array_minmax(sset_entity_count, Ioss::ParallelUtils::DO_MAX);
    }

    // If count > 0 on any processor, set the id if they have an id property.
    // If the containing sideset has an id and the sideblock doesn't have an id,
    // then set the sideblock id to the sideset id...
    size_t fb_index = 0;
    for (I=sidesets.begin(); I != sidesets.end(); ++I) {
      int fs_id = 0;
      if ((*I)->property_exists("id")) {
	fs_id = (*I)->get_property("id").get_int();
      }

      Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
      Ioss::SideBlockContainer::const_iterator J;

      if (fs_id > 0) {
	for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
	  if (!(*J)->property_exists("id")) {
	    (*J)->property_add(Ioss::Property("id", fs_id));
	  }
	  int count = sset_entity_count[fb_index++];
	  if (count > 0) {
	    set_id((*J), "Surface", &ids_);
	  }
	}
      } else {
	for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
	  int count = sset_entity_count[fb_index++];
	  if (count > 0) {
	    set_id((*J), "Surface", &ids_);
	  }
	}
      }
    }
    assert(fb_index == sset_entity_count.size());

    // The id has been set on all side blocks that had the id property.
    // Now, go through again and set the id on all side blocks.
    fb_index = 0;
    for (I=sidesets.begin(); I != sidesets.end(); ++I) {
      Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
      Ioss::SideBlockContainer::const_iterator J;

      for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
	int count = sset_entity_count[fb_index++];
	if (count > 0) {
	  get_id((*J), 'S', &ids_);
	  Ioxf::SideSet T(*(*J));
	  if (std::find(ssets.begin(), ssets.end(), T) == ssets.end()) {
	    ssets.push_back(T);
	  }
	} else {
	  // Set the "invalid" property.
	  Ioss::SideBlock *new_entity = const_cast<Ioss::SideBlock*>(*J);
	  new_entity->property_add(Ioss::Property("invalid", 1));
	}
      }
    }
    assert(fb_index == sset_entity_count.size());

    sidesetCount = ssets.size();

    // Nodesets ...
    {
      Ioss::NodeSetContainer nodesets = region->get_nodesets();
      Ioss::NodeSetContainer::const_iterator IN;
      for (IN=nodesets.begin(); IN != nodesets.end(); ++IN) {
	set_id(*IN, "Nodeset", &ids_);
      }

      for (IN=nodesets.begin(); IN != nodesets.end(); ++IN) {
	get_id(*IN, 'N', &ids_);
	const Ioxf::NodeSet T(*(*IN));
	if (std::find(nsets.begin(), nsets.end(), T) == nsets.end()) {
	  nsets.push_back(T);
	}
      }
      nodesetCount = nsets.size();
    }

    std::remove(hdfname.filename().c_str());  // does not matter if it fails

    InitXML(MainXML);

    Ioxf::CommunicationMetaData comm_meta;
    gather_communication_metadata(&comm_meta);

    {
      Ioss::SerializeIO	serializeIO__(this);

      *MainXML << "<BoundaryConditions>" << endl;
      WriteMetaXdmfNodesets(MainXML, nsets);
      WriteMetaXdmfSidesets(MainXML, ssets);
      *MainXML << "</BoundaryConditions>" << endl;
      WriteMetaXdmfElementBlock(blocks);

    }
  }

  void DatabaseIO::WriteMetaXdmfNodesets(std::ostringstream *XML, const std::vector<Ioxf::NodeSet> &nodesets)
  {
    int num_node_sets;
    num_node_sets = nodesets.size();

    if (num_node_sets){
      int node_set;
      int node_set_id;
      int num_nodes_in_set, num_df_in_set;

      for( node_set = 0 ; node_set < num_node_sets ; node_set++ ){
	node_set_id = nodesets[node_set].id;
	num_nodes_in_set = nodesets[node_set].nodeCount;
	num_df_in_set = nodesets[node_set].dfCount;
	NodeSetIndexNames.insert(BINMValuePair(std::string(nodesets[node_set].name),node_set));

	*XML << "\t<NodeSet Name=\"" << node_set_id << "\" >" << endl;
	*XML << "\t<Nodes>" << endl;
	*XML << "\t\t<DataStructure" << endl;
	*XML << "\t\t\tFormat=\"HDF\"" << endl;
	*XML << "\t\t\tDimensions=\"" << num_nodes_in_set << "\""<< endl;
	*XML << "\t\t\tDataType=\"Int\">" << endl;
	*XML << "\t\t\t" << hdfname.tailname()
	     << ":/NodeSets/" << node_set_id << "/Node" << endl;
	*XML << "\t\t</DataStructure>" << endl;
	*XML << "\t</Nodes>" << endl;
	if( num_df_in_set ){
	  *XML << "\t<DistributionFactor>" << endl;
	  *XML << "\t\t<DataStructure" << endl;
	  *XML << "\t\t\tFormat=\"HDF\"" << endl;
	  *XML << "\t\t\tDimensions=\"" << num_df_in_set << "\""<< endl;
	  *XML << "\t\t\tDataType=\"Float\">" << endl;
	  *XML << "\t\t\t" << hdfname.tailname()
	       << ":/NodeSets/" << node_set_id << "/DistributionFactor" << endl;
	  *XML << "\t\t</DataStructure>" << endl;
	  *XML << "\t</DistributionFactor>" << endl;
	}
	*XML << "\t</NodeSet>" << endl;
      }
    }
  }

  void DatabaseIO::WriteMetaXdmfSidesets(std::ostringstream *XML, const std::vector<Ioxf::SideSet> &sidesets)
  {
    int num_side_sets;
    num_side_sets = sidesets.size();
    if(num_side_sets){
      int side_set;
      int side_set_id;
      int num_sides_in_set, num_df_in_set;

      for( side_set = 0 ; side_set < num_side_sets ; side_set++ ){
	side_set_id = sidesets[side_set].id;
	num_sides_in_set = sidesets[side_set].sideCount;
	num_df_in_set = sidesets[side_set].dfCount;
	SideSetIndexNames.insert(BINMValuePair(std::string(sidesets[side_set].name),side_set));


	*XML << "\t<SideSet Name=\"" << side_set_id << "\" >" << endl;

	// there was no mention of this element block in the exodusII implementation
	// We might have to get rid of this.
	/*
	 *XML << "\t<Node>" << endl;
	 *XML << "\t\t<DataStructure" << endl;
	 *XML << "\t\t\tFormat=\"HDF\"" << endl;
	 *XML << "\t\t\tDimensions=\"" << SideSetNodes->GetNumberOfElements() << "\""<< endl;
	 *XML << "\t\t\tDataType=\"Int\">" << endl;
	 *XML << "\t\t\t" << hdfname.tailname()
	 << ":/SideSets/" << side_set_id << "/Node" << endl;
	 *XML << "\t\t</DataStructure>" << endl;
	 *XML << "\t</Node>" << endl;

	 num_element_in_set = sidesets[side_set].elemCount;
	 num_node_per_sideset= sidesets[side_set].nodesPerSideSet ;
	*/

	*XML << "\t<Side>" << endl;
	*XML << "\t\t<DataStructure" << endl;
	*XML << "\t\t\tFormat=\"HDF\"" << endl;
	*XML << "\t\t\tDimensions=\"" << num_sides_in_set << "\""<< endl;
	*XML << "\t\t\tDataType=\"Int\">" << endl;
	*XML << "\t\t\t" << hdfname.tailname()
	     << ":/SideSets/" << side_set_id << "/Side" << endl;
	*XML << "\t\t</DataStructure>" << endl;
	*XML << "\t</Side>" << endl;

	*XML << "\t<Element>" << endl;
	*XML << "\t\t<DataStructure" << endl;
	*XML << "\t\t\tFormat=\"HDF\"" << endl;
	*XML << "\t\t\tDimensions=\"" << num_sides_in_set << "\""<< endl;
	*XML << "\t\t\tDataType=\"Int\">" << endl;
	*XML << "\t\t\t" << hdfname.tailname()
	     << ":/SideSets/" << side_set_id << "/Element" << endl;
	*XML << "\t\t</DataStructure>" << endl;
	*XML << "\t</Element>" << endl;

	if( num_df_in_set ){
	  *XML << "\t<DistributionFactor>" << endl;
	  *XML << "\t\t<DataStructure" << endl;
	  *XML << "\t\t\tFormat=\"HDF\"" << endl;
	  *XML << "\t\t\tDimensions=\"" << num_df_in_set << "\""<< endl;
	  *XML << "\t\t\tDataType=\"Float\">" << endl;
	  *XML << "\t\t\t" << hdfname.tailname()
	       << ":/SideSets/" << side_set_id << "/DistributionFactor" << endl;
	  *XML << "\t\t</DataStructure>" << endl;
	  *XML << "\t</DistributionFactor>" << endl;
	}
	*XML << "\t</SideSet>" << endl;
      }
    }
  }
  void DatabaseIO::WriteMetaXdmfElementBlock(const std::vector<Ioxf::Block> &blocks)
  {

    char element_Type[MAX_STR_LENGTH+1];
    int num_elemets_blocks = blocks.size();

    std::ostringstream *XML = NULL;

    for ( int i = 0 ; i < num_elemets_blocks ; i++ ){
      int block_id = blocks[i].id;
      int num_element_in_block   = blocks[i].elementCount;
      int num_node_per_element = blocks[i].nodesPerElement;
      strcpy(element_Type,blocks[i].elType);

      // BlockIndexNames.insert(BINMValuePair(std::string(blocks[i].name),i));
      std::string block_name(blocks[i].name);
      BlockIndexNames.insert(BINMValuePair(block_name,i));

      XML = new std::ostringstream();
      BlockGridXmls.push_back(XML);


      *XML << "\t<Grid Name=\"" << block_name;
      if (isParallel) {
	*XML << " from Processor " << Ioss::Utils::to_string(myProcessor) << " of "
	     << Ioss::Utils::to_string(util().parallel_size()) << "\"";
	*XML << " Collection=\"" << block_name;
      }
      *XML << "\">" << endl;

      XML = new std::ostringstream();
      BlockParameterXmls.push_back(XML);

      // Parameter Iterations info will be defined later in end_state.

      XML = new std::ostringstream();
      BlockXmls.push_back(XML);


      *XML << "\t<Topology" << endl;
      if( strncasecmp( element_Type, "SHE", 3 ) == 0 ) {
	*XML << "\t\tType=\"Quadrilateral\""                              << endl;
      }
      else if( strncasecmp( element_Type, "HEX", 3 ) == 0 ) {
	*XML << "\t\tType=\"Hexahedron\""                                 << endl;
      }
      else if( strncasecmp( element_Type, "TET", 3 ) == 0 ) {
	*XML << "\t\tType=\"Tetrahedron\""                                << endl;
      }
      else if( strncasecmp( element_Type, "WED", 3 ) == 0 ) {
	*XML << "\t\tType=\"Wedge\""                                     << endl;
      }
      else if( strncasecmp( element_Type, "NOD", 3 ) == 0 ) {
	*XML << "\t\tType=\"PolyVertex\"" << endl;
      }
      else if( strncasecmp( element_Type, "BEA", 3 ) == 0 ) {
	*XML << "\t\tType=\"Polyline\""                                   << endl;
	*XML << "\t\tNodesPerElement=\"" <<  num_node_per_element << "\"" << endl;
      }
      else {
	*XML << "\t\tType=\"" << element_Type <<  "\""                    << endl;
      }
      *XML << "\t\tNumberOfElements=\"" << num_element_in_block << "\""         << endl;
      *XML << "\t\tBaseOffset=\"1\""                                            << endl;
      *XML << "\t\t>"                                                           << endl;
      if( strncasecmp( element_Type, "NOD", 3 ) != 0 ) {
	*XML << "\t\t<DataStructure"                                          << endl;
	*XML << "\t\t\tFormat=\"HDF\""                                        << endl;
	*XML << "\t\t\tDimensions=\"" << num_element_in_block << " "
	     << num_node_per_element << "\"" << endl;
	*XML << "\t\t\tDataType=\"Int\">"                                     << endl;
	*XML << "\t\t\t" << hdfname.tailname()
	     << ":/" << block_name << "/Connections"                          << endl;
	*XML << "\t\t</DataStructure>"                                        << endl;
      }
      *XML << "\t</Topology>"                                                   << endl;


      if (spatialDimension == 3) {
	*XML << "\t<Geometry Type=\"XYZ\">"                  << endl;
	*XML << "\t\t<DataStructure"                         << endl;
	*XML << "\t\t\tFormat=\"HDF\""                       << endl;
	*XML << "\t\t\tDimensions=\"" << nodeCount << " 3\"" << endl;
	*XML << "\t\t\tDataType=\"Float\">"                  << endl;
	*XML << "\t\t" << hdfname.tailname() << ":/Geometry" << endl;
	*XML << "\t\t</DataStructure>"                       << endl;
	*XML << "\t</Geometry>"                              << endl;
      }
      else {
	*XML << "\t<Geometry Type=\"XY\">"                   << endl;
	*XML << "\t\t<DataStructure"                         << endl;
	*XML << "\t\t\tFormat=\"HDF\""                       << endl;
	*XML << "\t\t\tDimensions=\"" << nodeCount << " 2\"" << endl;
	*XML << "\t\t\tDataType=\"Float\">"                  << endl;
	*XML << "\t\t" << hdfname.tailname() << ":/Geometry" << endl;
	*XML << "\t\t</DataStructure>"                       << endl;
	*XML << "\t</Geometry>"                              << endl;
      }

      XML = new std::ostringstream();
      BlockElementVarXmls.push_back(XML);
      // Element variables will be defined in  write_results_metadata

      XML = new std::ostringstream();
      BlockNodeVarXmls.push_back(XML);
      *XML << "\t<Attribute Name=\"" << "Values" << "\"" << endl;
      *XML << "\t\tCenter=\"" << "Node" << "\"" << endl;
      *XML << "\t\tType=\"" << "Scalar" << "\">" << endl;
      *XML << "\t\t<DataStructure" << endl;
      *XML << "\t\t\tFormat=\"HDF\"" << endl;
      *XML << "\t\t\tDimensions=\"" << nodeCount << "\"" << endl;
      *XML << "\t\t\tDataType=\""<< "Int" << "\" >" << endl;
      *XML << "\t\t" << hdfname.tailname() << ":/NodeData/ids" << endl;
      *XML << "\t\t</DataStructure>" << endl;
      *XML << "\t</Attribute>" << endl;

      XML = new std::ostringstream();
      BlockExtraAttributeXmls.push_back(XML);


      *XML << "\t<Attribute Name=\"BlockId\"" << endl;
      *XML << "\t\tCenter=\"Grid\""               << endl;
      *XML << "\t\tType=\"Scalar\">"              << endl;
      *XML << "\t\t<DataStructure"                << endl;
      *XML << "\t\t\tFormat=\"XML\""              << endl;
      *XML << "\t\t\tDimensions=\"1\""            << endl;
      *XML << "\t\t\tDataType=\"Int\" >"          << endl;
      *XML << "\t\t" << block_id         << endl;
      *XML << "\t\t</DataStructure>"              << endl;
      *XML << "\t</Attribute>"                    << endl;


      XML = new std::ostringstream();
      BlockFinishXmls.push_back(XML);
      // Grid & Parameter Close statemets will be defined in end_state.

    }
  }

  void DatabaseIO::WriteXmlFile(int nIterations)
  {

    std::ofstream XMLFile(xmlname.filename().c_str());
    std::ostringstream *XML = NULL;

    //  MainXML = FixXML(MainXML);

    XMLFile << MainXML->str() << endl;
    if (!XMLFile.good() || XMLFile.fail() || XMLFile.bad())
      xdmf_error(xmlname.filename().c_str(), __LINE__ , myProcessor);


    int num_elemets_blocks;
    num_elemets_blocks = BlockGridXmls.size();
    for (int i = 0 ; i < num_elemets_blocks ; i++ ){

      XML = BlockGridXmls[i];
      XMLFile << XML->str() << endl;

      XML = BlockParameterXmls[i];
      XMLFile << XML->str() << endl;

      if( nIterations > 0 ){
	XMLFile << "\t<Parameter Name=\"Iteration\" Type=\"Range\"" << endl;
	XMLFile << "\t\tFormat=\"Iteration %d\" Values=\"" << 1 <<
	  " 1 " << nIterations << "\"" << endl;
	XMLFile << "\t\tCurrentIndex=\"0\">" << endl;
      }


      XML = BlockXmls[i];
      XMLFile << XML->str() << endl;

      XML = BlockElementVarXmls[i];
      XMLFile << XML->str() << endl;


      XML = BlockNodeVarXmls[i];
      XMLFile << XML->str() << endl;

      XML = BlockExtraAttributeXmls[i];
      XMLFile << XML->str() << endl;

      XML = BlockFinishXmls[i];
      XMLFile << XML->str() << endl;
      if( nIterations > 0 ){
	XMLFile << "\t</Parameter>" << endl;
      }
      XMLFile << "</Grid>" << endl;
      XMLFile << "" << endl;

    }
    XMLFile << " </Domain>" << endl<< "</Xdmf>" << endl;
    // XMLFile << '\0';

    if (!XMLFile.good() || XMLFile.fail() || XMLFile.bad())
      xdmf_error(xmlname.filename().c_str(), __LINE__ , myProcessor);

    XMLFile.close();



    Ioss::SerializeIO	serializeIO__(this);

    MergeXmlFiles();
  }

  std::string DatabaseIO::decode_proc_filename(const std::string &filename,  int processor)
  {
    std::string decoded_filename(filename);

    if (isParallel) {
      // Running in parallel, assume nemesis and decode what the filename
      // should be for this processor.
      int num_processors = util().parallel_size();

      // Current format for per-processor file names is:
      // PREFIX/basename.num_proc.cur_proc
      // the 'cur_proc' field is padded to be the same width as
      // the 'num_proc' field
      // Examples: basename.8.1, basename.64.03, basename.128.001

      // Create a std::string containing the total number of processors
      std::string num_proc = Ioss::Utils::to_string(num_processors);
      int proc_width = num_proc.length();

      // Create a std::string containing the current processor number
      std::string cur_proc = Ioss::Utils::to_string(processor);
      int cur_width = cur_proc.length();

      // Build the filename
      // *FIX* Ignore raid issues for now
      decoded_filename += ".";
      decoded_filename += num_proc;
      decoded_filename += ".";

      // Now, pad with zeros so that 'cur_proc' portion is same
      // width as 'num_proc' portion.
      while (cur_width++ < proc_width) {
	decoded_filename += "0";
      }

      decoded_filename += cur_proc;
    }
    return decoded_filename;
  }


  void DatabaseIO::MergeXmlFiles()
  {
    Ioss::FileInfo filename(get_filename()+".xmf");
    // std::ostringstream *XML = NULL;

#ifdef HAVE_MPI
    MPI_Barrier(util().communicator());
#endif
    if ((myProcessor == 0) && isParallel) {

      xmlDocPtr newdoc = xmlNewDoc((xmlChar*)"1.0");
      xmlNodePtr newxdmfroot = newdoc->children = xmlNewDocNode(newdoc, NULL, (xmlChar*) "Xdmf", NULL);
      xmlNodePtr newxdmfdomain = xmlNewChild(newxdmfroot, NULL, (xmlChar*)"Domain", NULL);

      // xmlSetProp(newroot, (xmlChar*)"version", (xmlChar*)"1.0");

      for (int i = 0 ; i <  util().parallel_size() ; i++ ){
	std::string docname = decode_proc_filename(std::string(filename.basename()),i) + ".xmf";
	xmlDocPtr doc = xmlParseFile(docname.c_str());
	if (doc == NULL)
	  xdmf_error(docname, __LINE__ , myProcessor);

	xmlNodePtr cur = xmlDocGetRootElement(doc);
	if (cur == NULL)
	  xdmf_error(docname + " bad xml file", __LINE__ , myProcessor);

	if (xmlStrcmp(cur->name, (const xmlChar*) "Xdmf"))
	  xdmf_error(docname + " missing Xdmf block", __LINE__ , myProcessor);

	cur = cur->xmlChildrenNode;
	while ((cur != NULL) &&  (xmlStrcmp(cur->name, (const xmlChar*) "Domain")) )
	  cur = cur->next;

	cur = cur->xmlChildrenNode;
	while( cur != NULL) {
	  if (!xmlStrcmp(cur->name, (const xmlChar*) "Grid")) {
	    xmlNodePtr node_copy = xmlCopyNode(cur, 1);
	    xmlAddChild(newxdmfdomain, node_copy); // append grid to new document
	  }
	  cur =  cur->next ;
	}
	xmlFreeDoc(doc);

      }
      xmlIndentTreeOutput = 1;
      int ret_val = xmlSaveFormatFileEnc(filename.filename().c_str(), newdoc, NULL, 1);
      if (ret_val == -1)
	xdmf_error(" Unable to create " + std::string(filename.filename()), __LINE__ , myProcessor);

      xmlFreeDoc(newdoc);

    }
  }


  void DatabaseIO::gather_communication_metadata(Ioxf::CommunicationMetaData *meta)
  {
    if (isParallel) {

      meta->processorCount = util().parallel_size();
      meta->processorId = myProcessor;

      if (get_region()->property_exists("global_node_count")) {
	meta->globalNodes = get_region()->
	  get_property("global_node_count").get_int();
      }

      if (get_region()->property_exists("global_element_count")) {
	meta->globalElements = get_region()->
	  get_property("global_element_count").get_int();
      }

      // ========================================================================
      // Load balance parameters (NEMESIS, p15)
      meta->nodesInternal = nodeCount;
      meta->nodesBorder   = 0;
      meta->nodesExternal = 0; // Shadow nodes == 0 for now
      meta->elementsInternal = elementCount;
      meta->elementsBorder   = 0;

      // Now, see if any of the above are redefined by a property...
      if (get_region()->property_exists("internal_node_count")) {
	meta->nodesInternal = get_region()->
	  get_property("internal_node_count").get_int();
      }

      if (get_region()->property_exists("border_node_count")) {
	meta->nodesBorder = get_region()->
	  get_property("border_node_count").get_int();
      }

      if (get_region()->property_exists("internal_element_count")) {
	meta->elementsInternal = get_region()->
	  get_property("internal_element_count").get_int();
      }

      if (get_region()->property_exists("border_element_count")) {
	meta->elementsBorder = get_region()->
	  get_property("border_element_count").get_int();
      }

      Ioss::CommSetContainer comm_sets = get_region()->get_commsets();
      Ioss::CommSetContainer::const_iterator I = comm_sets.begin();
      while (I != comm_sets.end()) {

	Ioss::CommSet *cs = *I;

	std::string type = cs->get_property("entity_type").get_string();
	int count = cs->get_property("entity_count").get_int();
	int id = get_id(cs, 'C', &ids_);

	if (type == "node") {
	  meta->nodeMap.push_back(Ioxf::CommunicationMap(id, count, 'n'));
	} else if (type == "side") {
	  meta->elementMap.push_back(Ioxf::CommunicationMap(id, count, 'e'));
	} else {
	  std::ostringstream errmsg;
	  errmsg << "INTERNAL ERROR in gather_communication_data";
	  IOSS_ERROR(errmsg);
	}
	++I;
      }
    }
    commsetNodeCount = meta->nodeMap.size();
    commsetElemCount = meta->elementMap.size();
  }

  void DatabaseIO::add_region_fields() const
  {
    int field_count = add_results_fields("G", get_region(), 1);
    DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
    new_this->globalValues.resize(field_count);
  }

  int DatabaseIO::add_results_fields(char const */* type */,
				     Ioss::GroupingEntity */* entity */,
				     int /* count */, int /* position */) const
  {
    output_only();
    return 0;
  }

  Ioss::Field DatabaseIO::get_next_field(char** /* names */, int */* index */,
					 int /* num_names */, int /* count */,
					 int */* truth_table */) const
  {
    return Ioss::Field("", Ioss::Field::REAL, "scalar",
		       Ioss::Field::TRANSIENT, 0);
  }

  void DatabaseIO::write_results_metadata()
  {
    // Does not yet support results on sideblocks or any sets
    int glob_index = 0;
    glob_index = gather_names("G", get_region(), glob_index, true);
    assert(glob_index == (int)globalVariables.size());

    {
      Ioss::NodeBlockContainer node_blocks = get_region()->get_node_blocks();
      assert(node_blocks.size() == 1);
      Ioss::NodeBlockContainer::const_iterator I;

      int index = 0;
      for (I=node_blocks.begin(); I != node_blocks.end(); ++I) {
	glob_index = gather_names("N", *I, glob_index, true);
	index = gather_names("N", *I, index, false);
      }
      assert(index == (int)nodalVariables.size());
    }

    {
      Ioss::ElementBlockContainer element_blocks =
	get_region()->get_element_blocks();
      assert(check_block_order(element_blocks));
      Ioss::ElementBlockContainer::const_iterator I;

      int index = 0;
      for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
	glob_index = gather_names("E", *I, glob_index, true);
	index = gather_names("E", *I, index, false);
	// Erik example of how to pull out the name
	// std::string ge_name = (*I)->name();
      }
      assert(index == (int)elementVariables.size());
      generate_var_xmltable("E");
    }
    //    Not doing anything with Nodeset and sideset Variables yet....

    {
      Ioss::NodeSetContainer nodesets = get_region()->get_nodesets();
      Ioss::NodeSetContainer::const_iterator I;

      int index = 0;
      for (I=nodesets.begin(); I != nodesets.end(); ++I) {
	glob_index = gather_names("M", *I, glob_index, true);
	index = gather_names("M", *I, index, false);
      }
      assert(index == (int)nodesetVariables.size());
      generate_var_xmltable("M");
    }

    {
      int index = 0;
      Ioss::SideSetContainer sidesets = get_region()->get_sidesets();
      Ioss::SideSetContainer::const_iterator I;
      for (I=sidesets.begin(); I != sidesets.end(); ++I) {
	Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
	Ioss::SideBlockContainer::const_iterator J;

	for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
	  glob_index = gather_names("S", *J, glob_index, true);
	  index = gather_names("S", *J, index, false);
	}
      }
      assert(index == (int)sidesetVariables.size());
      generate_var_xmltable("S");
    }
  }

  int DatabaseIO::gather_names(const char *type,
			       const Ioss::GroupingEntity *ge,
			       int index, bool reduction)
  {
    DatabaseIO *new_this = const_cast<DatabaseIO*>(this);
    int new_index = index;
    assert(type[0] == 'E' || type[0] == 'N' || type[0] == 'G' || type[0] == 'M' || type[0] == 'S');

    bool nblock = (type[0] == 'N');

    // Get names of all transient and reduction fields...
    Ioss::NameList results_fields;
    if (reduction)
      ge->field_describe(Ioss::Field::REDUCTION, &results_fields);

    if (!reduction || type[0] == 'G')
      ge->field_describe(Ioss::Field::TRANSIENT, &results_fields);

    // NOTE: For exodusII, the convention is that the displacement
    //       fields are the first 'ndim' fields in the file.
    //       Try to find a likely displacement field
    std::string disp_name;
    bool has_disp = false;
    if (!reduction && nblock && new_index == 0) {
      has_disp = find_displacement_field(results_fields, ge, spatialDimension,
					 &disp_name);
      if (has_disp)
	new_index += spatialDimension;
    }

    int save_index = 0;
    Ioss::NameList::const_iterator IF;

    for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
      std::string field_name = *IF;

      if (has_disp && field_name == disp_name && new_index != 0) {
	save_index = new_index;
	new_index = 0;
      }

      Ioss::Field field = ge->get_field(field_name);
      const Ioss::VariableType *var_type = field.transformed_storage();

      for (int i=1; i <= var_type->component_count(); i++) {
	std::string var_string = var_type->label_name(field_name, i, get_field_separator());

	// Add to 'VariableNameMap elementVariables' or
	// 'VariableNameMap nodalVariables' so can determine exodusII
	// index given a Sierra field name.  exodusII index is just 'i+1'
	if (reduction || type[0] == 'G') {
	  // If this is not a global (region) variable, need to prepend the block name
	  // to avoid name collisions...
	  if (type[0] == 'E' || type[0] == 'M' || type[0] == 'S') {
	    std::string ge_name = ge->name();
	    var_string = ge_name + ":" + var_string;
	  }
	  if (globalVariables.find(var_string) == globalVariables.end()) {
	    new_this->globalVariables.insert(VNMValuePair(var_string,
							  ++new_index));
	  }
	} else if (type[0] == 'N') {
	  if (nodalVariables.find(var_string) == nodalVariables.end()) {
	    new_this->nodalVariables.insert(VNMValuePair(var_string,
							 ++new_index));
	  }
	} else if (type[0] == 'E') {
	  if (elementVariables.find(var_string) == elementVariables.end()) {
	    new_this->elementVariables.insert(VNMValuePair(var_string,
							   ++new_index));
	  }
	} else if (type[0] == 'M') {
	  if (nodesetVariables.find(var_string) == nodesetVariables.end()) {
	    new_this->nodesetVariables.insert(VNMValuePair(var_string,
							   ++new_index));
	  }
	} else if (type[0] == 'S') {
	  if (sidesetVariables.find(var_string) == sidesetVariables.end()) {
	    new_this->sidesetVariables.insert(VNMValuePair(var_string,
							   ++new_index));
	  }
	}
      }
      if (has_disp && field_name == disp_name) {
	new_index = save_index;
      }
    }
    return new_index;
  }

  void DatabaseIO::generate_var_xmltable(const char *type)
  {
    if (type[0] == 'E')
      generate_element_var_xmltable();
    else if (type[0] == 'M')
      generate_nodeset_var_xmltable();
    else if (type[0] == 'S')
      generate_sideset_var_xmltable();
  }

  int DatabaseIO::get_xml_stream(const std::string &block_name) const
  {
    int BlockIndex = -1;
    IndexNameMap::iterator BN = BlockIndexNames.find(block_name);
    
    if (BN != BlockIndexNames.end()) {
      BlockIndex =  (*BN).second;
    }
    
    return BlockIndex ;
  }

  void DatabaseIO::generate_element_var_xmltable()
  {
    int var_count = elementVariables.size();

    if (var_count == 0 || elementBlockCount == 0)
      return;

    // Member variable.  Will be deleted in destructor...
    elementTruthTable = new int[elementBlockCount*var_count];

    // Zero-fill truth table; set to '1' when variable found...
    std::fill_n(&elementTruthTable[0], elementBlockCount*var_count, 0);

    // Fill in the truth table.  It is conceptually a two-dimensional array of
    // the form 'array[num_element_blocks][num_element_var]'.  In C++,
    // the values for the first element block are first, followed by
    // next element block, ...
    Ioss::ElementBlockContainer element_blocks = get_region()->get_element_blocks();
    Ioss::ElementBlockContainer::const_iterator I;

    assert(check_block_order(element_blocks));

    int offset = 0;
    for (I=element_blocks.begin(); I != element_blocks.end(); ++I) {
      // Get names of all transient and reduction fields...
      Ioss::NameList results_fields;
      (*I)->field_describe(Ioss::Field::TRANSIENT, &results_fields);
      (*I)->field_describe(Ioss::Field::REDUCTION, &results_fields);

      std::string block_name = (*I)->name();
      int count = (*I)->get_property("entity_count").get_int();

      std::ostringstream *XML=NULL;
      int BlockIndex = get_xml_stream(block_name);
      XML = BlockElementVarXmls[BlockIndex];

      Ioss::NameList::const_iterator IF;
      for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
	std::string field_name = *IF;

	Ioss::Field field = (*I)->get_field(field_name);
	const Ioss::VariableType *var_type = field.transformed_storage();

	for (int i=1; i <= var_type->component_count(); i++) {
	  std::string var_string = var_type->label_name(field_name, i, get_field_separator());

	  // Find position of 'var_string' in 'elementVariables'
	  VariableNameMap::iterator VN = elementVariables.find(var_string);
	  if (XML && VN != elementVariables.end()) {
	    *XML << "\t<Attribute Name=\"" << var_string << "\""  << endl;
	    *XML << "\t\tCenter=\"" << "Cell" << "\""             << endl;
	    *XML << "\t\tType=\"" << "Scalar" << "\">"            << endl;
	    *XML << "\t\t<DataStructure"                          << endl;
	    *XML << "\t\t\tFormat=\"HDF\""                        << endl;
	    *XML << "\t\t\tDimensions=\"" << count << "\""        << endl;
	    *XML << "\t\t\tDataType=\""<< "Float"  << "\">"       << endl;
	    *XML << "\t\t" << hdfname.tailname()
		 << ":/" << block_name << "/Iteration/" << var_string << endl;
	    *XML << "\t\t</DataStructure>"                        << endl;
	    *XML << "\t</Attribute>"                              << endl;
	  }
	}
      }
      offset += var_count;
    }
    assert(offset == var_count * elementBlockCount);
  }

  void DatabaseIO::generate_sideset_var_xmltable()
  {
    int var_count = sidesetVariables.size();

    if (var_count == 0 || sidesetCount == 0)
      return;

    // Member variable.  Will be deleted in destructor...
    sidesetTruthTable = new int[sidesetCount*var_count];

    // Zero-fill truth table; set to '1' when variable found...
    std::fill_n(&sidesetTruthTable[0], sidesetCount*var_count, 0);

    // Fill in the truth table.  It is conceptually a two-dimensional array of
    // the form 'array[num_blocks][num_var]'.  In C++,
    // the values for the first block are first, followed by
    // next block, ...
    int offset = 0;

    Ioss::SideSetContainer sidesets = get_region()->get_sidesets();
    Ioss::SideSetContainer::const_iterator I;

    for (I=sidesets.begin(); I != sidesets.end(); ++I) {
      Ioss::SideBlockContainer side_blocks = (*I)->get_side_blocks();
      Ioss::SideBlockContainer::const_iterator J;

      bool found_one = false;
      for (J=side_blocks.begin(); J != side_blocks.end(); ++J) {
	// See if this sideblock has a corresponding entry in the sideset list.
	if ((*J)->property_exists("invalid"))
	  continue;

	found_one = true;
	// Get names of all transient and reduction fields...
	Ioss::NameList results_fields;
	(*J)->field_describe(Ioss::Field::TRANSIENT, &results_fields);
	(*J)->field_describe(Ioss::Field::REDUCTION, &results_fields);

	Ioss::NameList::const_iterator IF;
	for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
	  std::string field_name = *IF;

	  Ioss::Field field = (*J)->get_field(field_name);
	  const Ioss::VariableType *var_type = field.transformed_storage();

	  for (int i=1; i <= var_type->component_count(); i++) {
	    std::string var_string = var_type->label_name(field_name, i, get_field_separator());
	    // Find position of 'var_string' in 'elementVariables'
	    VariableNameMap::iterator VN = sidesetVariables.find(var_string);
	    if (VN != sidesetVariables.end()) {
	      // Index '(*VN).second' is 1-based...
	      sidesetTruthTable[offset + (*VN).second-1] = 1;
	    }
	  }
	}
      }
      if (found_one)
	offset += var_count;
    }
  }

  void DatabaseIO::generate_nodeset_var_xmltable()
  {
    int var_count = nodesetVariables.size();

    if (var_count == 0 || nodesetCount == 0)
      return;

    // Member variable.  Will be deleted in destructor...
    nodesetTruthTable = new int[nodesetCount*var_count];

    // Zero-fill truth table; set to '1' when variable found...
    std::fill_n(&nodesetTruthTable[0], nodesetCount*var_count, 0);

    // Fill in the truth table.  It is conceptually a two-dimensional array of
    // the form 'array[num_blocks][num_var]'.  In C++,
    // the values for the first block are first, followed by
    // next block, ...
    Ioss::NodeSetContainer nodesets = get_region()->get_nodesets();
    Ioss::NodeSetContainer::const_iterator I;

    int offset = 0;
    for (I=nodesets.begin(); I != nodesets.end(); ++I) {
      // Get names of all transient and reduction fields...
      Ioss::NameList results_fields;
      (*I)->field_describe(Ioss::Field::TRANSIENT, &results_fields);
      (*I)->field_describe(Ioss::Field::REDUCTION, &results_fields);

      Ioss::NameList::const_iterator IF;
      for (IF = results_fields.begin(); IF != results_fields.end(); ++IF) {
	std::string field_name = *IF;

	Ioss::Field field = (*I)->get_field(field_name);
	const Ioss::VariableType *var_type = field.transformed_storage();

	for (int i=1; i <= var_type->component_count(); i++) {
	  std::string var_string = var_type->label_name(field_name, i, get_field_separator());
	  // Find position of 'var_string' in 'elementVariables'
	  VariableNameMap::iterator VN = nodesetVariables.find(var_string);
	  if (VN != nodesetVariables.end()) {
	    // Index '(*VN).second' is 1-based...
	    nodesetTruthTable[offset + (*VN).second-1] = 1;
	  }
	}
      }
      offset += var_count;
    }
    assert(offset == var_count * nodesetCount);
  }

  void DatabaseIO::build_element_reorder_map(int start, int count)
  {
    // Note: To further add confusion, the reorderElementMap is 0-based
    // and the reverseElementMap and elementMap are 1-based. This is
    // just a consequence of how they are intended to be used...
    //
    // start is based on a 0-based array -- start of the reorderMap to build.

    if (reorderElementMap.empty())
      reorderElementMap.resize(elementCount);
    assert(static_cast<int>(elementMap.size()) == elementCount+1);
    assert(static_cast<int>(reverseElementMap.size()) <= elementCount); // Built in pieces

    int end = start+count;
    for (int i=start; i < end; i++) {
      int global_id = elementMap[i+1];
      int orig_local_id = element_global_to_local(global_id) - 1;

      // If we assume that partial output is not being used (it
      // currently isn't in Sierra), then the reordering should only be
      // a permutation of the original ordering within this element block...
      assert(orig_local_id >= start && orig_local_id <= end);
      reorderElementMap[i] = orig_local_id;
    }
  }

  void DatabaseIO::build_node_reorder_map(int *new_ids, int count)
  {
    // This routine builds a map that relates the current node id order
    // to the original node ordering in affect at the time the file was
    // created. That is, the node map used to define the topology of the
    // model.  Now, if there are changes in node ordering at the
    // application level, we build the node reorder map to map the
    // current order into the original order.  An added complication is
    // that this is more than just a reordering... It may be that the
    // application has 'ghosted' nodes that it doesnt want put out on
    // the database, so the reorder map must handle a node that is not
    // in the original mesh and map that to an invalid value (currently
    // using -1 as invalid value...)


    // Note: To further add confusion,
    // the reorderNodeMap and new_ids are 0-based
    // the reverseNodeMap and nodeMap are 1-based. This is
    // just a consequence of how they are intended to be used...

    reorderNodeMap.resize(count);

    for (int i=0; i < count; i++) {
      int global_id = new_ids[i];

      // This will return 0 if node is not found in list.
      int orig_local_id = node_global_to_local(global_id, false) - 1;

      reorderNodeMap[i] = orig_local_id;
    }
  }

  // Handle special output time requests -- primarily restart (cycle, keep, overwrite)
  // Given the global region step, return the step on the database...
  int DatabaseIO::get_database_step(int global_step) const
  {
    assert(overlayCount >= 0 && cycleCount >= 0);
    if (overlayCount == 0 && cycleCount == 0)
      return global_step;

    int local_step = global_step - 1;
    local_step /= (overlayCount + 1);
    if (cycleCount > 0)
      local_step %= cycleCount;
    return local_step + 1;
  }
}

namespace {
bool set_id(const Ioss::GroupingEntity *entity, const char *type, Ioxf::EntityIdSet *idset)
{
  // See description of 'get_id' function.  This function just primes
  // the idset with existing ids so that when we start generating ids,
  // we don't overwrite an existing one.

  // Avoid a few string constructors/destructors
  static std::string prop_name("name");
  static std::string id_prop("id");

  bool succeed = false;
  if (entity->property_exists(id_prop)) {
    int id = entity->get_property(id_prop).get_int();

    // See whether it already exists...
    succeed = idset->insert(std::make_pair(type[0],id)).second;
    if (!succeed) {
      std::string name_string = entity->get_property(prop_name).get_string();
      IOSS_WARNING << "The " << type << " named '" << name_string
		       << "' has the exodus id = " << id
		       << "\nwhich has already been assigned to another " << type
		       << ".  This entity will be assigned a new id.\n"
		       << "Contact gdsjaar@sandia.gov if you need explanation of this warning.";

      // Need to remove the property so it doesn't cause problems
      // later...
      Ioss::GroupingEntity *new_entity = const_cast<Ioss::GroupingEntity*>(entity);
      new_entity->property_erase(id_prop);
      assert(!entity->property_exists(id_prop));
    }
  }
  return succeed;
}

int get_id(const Ioss::GroupingEntity *entity, char type, Ioxf::EntityIdSet *idset)
{
  // Sierra uses names to refer to grouping entities; however,
  // exodusII requires integer ids.  When reading an exodusII file,
  // the DatabaseIO creates a name by concatenating the entity
  // type (e.g., 'block') and the id separated by an underscore.  For
  // example, an exodusII element block with an id of 100 would be
  // encoded into "block_100"

  // This routine tries to determine the id of the entity using 3
  // approaches:
  //
  // 1. If the entity contains a property named 'id', this is used.
  // The DatabaseIO actually stores the id in the "id" property;
  // however, other grouping entity creators are not required to do
  // this so the property is not guaranteed to exist.
  //
  // 2.If property does not exist, it tries to decode the entity name
  // based on the above encoding.  Again, it is not required that the
  // name follow this convention so success is not guaranteed.
  //
  // 3. If all other schemes fail, the routine picks an id for the entity
  // and returns it.  It also stores this id in the "id" property so an
  // entity will always return the same id for multiple calls.
  // Note that this violates the 'const'ness of the entity so we use
  // a const-cast.

  // Avoid a few string constructors/destructors
  static std::string prop_name("name");
  static std::string id_prop("id");

  int id = 1;

  if (entity->property_exists(id_prop)) {
    id = entity->get_property(id_prop).get_int();
    return id;

  } else {
    // Try to decode an id from the name.
    std::string name_string = entity->get_property(prop_name).get_string();
    const char *name = name_string.c_str();
    int len = std::strlen(name);

    // Start 1 character from end just to avoid problems with a trailing '_'
    // Assume that if '_' found, there is at least 1 character following.
    for (int i=len-2; i >= 0; i--) {
      if (name[i] == '_') {
	id = std::atoi(&name[i+1]);
	break;
      } else if (!isdigit(name[i])) {
	break;
      }
    }
  }

  // At this point, we either have an id equal to '1' or we have an id
  // extracted from the entities name. Increment it until it is
  // unique...
  while (idset->find(std::make_pair(type, id)) != idset->end()) {
    ++id;
  }

  // 'id' is a unique id for this entity type...
  idset->insert(std::make_pair(type,id));
  Ioss::GroupingEntity *new_entity = const_cast<Ioss::GroupingEntity*>(entity);
  new_entity->property_add(Ioss::Property(id_prop, id));
  return id;
}

bool find_displacement_field(Ioss::NameList &fields,
			     const Ioss::GroupingEntity *block,
			     int ndim,
			     std::string *disp_name)
{
  // This is a kluge to work with many of the SEACAS codes.  The
  // convention used (in Blot and others) is that the first 'ndim'
  // nodal variables are assumed to be displacements *if* the first
  // character of the names is 'D' and the last characters match the
  // coordinate labels (typically 'X', 'Y', and 'Z').  This routine
  // looks for the first field that begins with 'd' and is of the
  // correct storage type (VECTOR_2D or VECTOR_3D).  If found, it
  // returns the name.
  //
  // Sometime this should check all names and return the 'most
  // likely...'

  Ioss::NameList::const_iterator IF;
  Ioss::NameList::const_iterator IFend = fields.end();
  for (IF = fields.begin(); IF != IFend; ++IF) {
    const char *name = (*IF).c_str();

    if (name[0] == 'd' || name[0] == 'D') {
      const Ioss::VariableType *var_type =
	block->get_field((*IF)).transformed_storage();
      int comp_count = var_type->component_count();
      if (comp_count == ndim) {
	*disp_name = *IF;
	return true;
      }
    }
  }
  return false;
}

int field_warning(const Ioss::GroupingEntity *ge,
		  const Ioss::Field &field, const std::string& inout)
{
  IOSS_WARNING << ge->type() << " '"
		   << ge->name()
		   << "'. Unknown " << inout << " field '"
		   << field.get_name() << "'";
  return -4;
}

#ifndef NDEBUG
bool check_block_order(const Ioss::ElementBlockContainer &blocks)
{
  // Verify that element blocks are defined in sorted offset order...
  Ioss::ElementBlockContainer::const_iterator I;

  int eb_offset = -1;
  for (I=blocks.begin(); I != blocks.end(); ++I) {
    int this_off = (*I)->get_offset();
    if (this_off < eb_offset)
      return false;
    eb_offset = this_off;
  }
  return true;
}
#endif
}
