// CGNS Assumptions:
// * All boundary conditions are listed as Family nodes at the "top" level.
// * Unstructured mesh only
// * Single element block per zone.
// * Serial for now.
// * Single Base.
// * ZoneGridConnectivity is 1to1 with point lists

// Copyright(C) 2015
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
#include <Ioss_Utils.h>
#include <assert.h>
#include <cgns/Iocgns_DatabaseIO.h>
#include <stddef.h>
#include <sys/select.h>
#include <time.h>
#include <numeric>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cgnslib.h>

#if !defined(CGNSLIB_H)
#error "Could not include cgnslib.h"
#endif

#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_EntityType.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"

#include "Ioss_Field.h"
#include "Ioss_IOFactory.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_State.h"
#include "Ioss_Utils.h"
#include "Ioss_VariableType.h"

namespace {
  const char *Version() {return "Iocgns_DatabaseIO.C 2016/01/28";}

  void cgns_error(int cgnsid, int lineno, int /* processor */)
  {
    std::ostringstream errmsg;
    errmsg << "CGNS error '" << cg_get_error() << "' at line " << lineno
	   << " in file '" << Version()
	   << "' Please report to gdsjaar@sandia.gov if you need help.";
    if (cgnsid > 0) {
      cg_close(cgnsid);
    }
    IOSS_ERROR(errmsg);
  }

  template <typename INT>
  void map_cgns_face_to_ioss(const Ioss::ElementTopology *parent_topo, size_t num_to_get, INT *idata)
  {
    // The {topo}_map[] arrays map from CGNS face# to IOSS face#.
    // See http://cgns.github.io/CGNS_docs_current/sids/conv.html#unstructgrid
    // NOTE: '0' for first entry is to account for 1-based face numbering.

    switch (parent_topo->shape())
      {
      case Ioss::ElementShape::HEX:
	static int hex_map[] = {0, 5, 1, 2, 3, 4, 6};
	for (size_t i=0; i < num_to_get; i++) {
	  idata[2*i+1] = hex_map[idata[2*i+1]];
	}
	break;

      case Ioss::ElementShape::TET:
	static int tet_map[] = {0, 4, 1, 2, 3};
	for (size_t i=0; i < num_to_get; i++) {
	  idata[2*i+1] = tet_map[idata[2*i+1]];
	}
	break;

      case Ioss::ElementShape::PYRAMID:
	static int pyr_map[] = {0, 5, 1, 2, 3, 4};
	for (size_t i=0; i < num_to_get; i++) {
	  idata[2*i+1] = pyr_map[idata[2*i+1]];
	}
	break;

      case Ioss::ElementShape::WEDGE:
#if 0
	static int wed_map[] = {0, 1, 2, 3, 4, 5}; // Same
	// Not needed -- maps 1 to 1
	for (size_t i=0; i < num_to_get; i++) {
	  idata[2*i+1] = wed_map[idata[2*i+1]];
	}
#endif
	break;
      default:
	;
      }
  }

  std::string map_cgns_to_topology_type(CG_ElementType_t type)
  {
    std::string topology = "unknown";
    switch (type)
      {
      case CG_NODE:
	topology = "tetra4"; break;
      case CG_BAR_2:
	topology = "bar2"; break;
      case CG_BAR_3:
	topology = "bar3"; break;
      case CG_TRI_3:
	topology = "tri3"; break;
      case CG_TRI_6:
	topology = "tri6"; break;
      case CG_QUAD_4:
	topology = "quad4"; break;
      case CG_QUAD_8:
	topology = "quad8"; break;
      case CG_QUAD_9:
	topology = "quad9"; break;
      case CG_TETRA_4:
	topology = "tetra4"; break;
      case CG_TETRA_10:
	topology = "tetra10"; break;
      case CG_PYRA_5:
	topology = "pyramid5"; break;
      case CG_PYRA_13:
	topology = "pyramid13"; break;
      case CG_PYRA_14:
	topology = "pyramid14"; break;
      case CG_PENTA_6:
	topology = "wedge6"; break;
      case CG_PENTA_15:
	topology = "wedge15"; break;
      case CG_PENTA_18:
	topology = "wedge18"; break;
      case CG_HEXA_8:
	topology = "hex8"; break;
      case CG_HEXA_20:
	topology = "hex20"; break;
      case CG_HEXA_27:
	topology = "hex27"; break;
      default:
	std::cerr << "WARNING: Found topology of type "
		  << cg_ElementTypeName(type)
		  << " which is not currently supported.\n";
	topology = "unknown";
      }
    return topology;
  }
}
namespace Iocgns {

  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
			 Ioss::DatabaseUsage db_usage,
			 MPI_Comm communicator,
			 const Ioss::PropertyManager &props) :
    Ioss::DatabaseIO(region, filename, db_usage, communicator, props),
    cgnsFilePtr(-1)
  {
    dbState = Ioss::STATE_UNKNOWN;

    if (!is_input()) {
      std::ostringstream errmsg;
      errmsg << "ERROR: IOSS CGNS Currently only supports reading, not writing.\n";
      IOSS_ERROR(errmsg);
    }

    std::cout << "CGNS DatabaseIO using " << CG_SIZEOF_SIZE << "-bit integers.\n";
    if (CG_SIZEOF_SIZE == 64) {
      set_int_byte_size_api(Ioss::USE_INT64_API);
    }

    openDatabase();
  }

  void DatabaseIO::openDatabase() const
  {
    if (cgnsFilePtr < 0) {
      if (is_input()) {
        int ierr = cg_open(get_filename().c_str(), CG_MODE_READ, &cgnsFilePtr);
	if (ierr != CG_OK) {
	  // NOTE: Code will not continue past this call...
	  std::ostringstream errmsg;
	  errmsg << "ERROR: Problem opening file '" << get_filename() << "' for read access.";
	  IOSS_ERROR(errmsg);
	}
      }
    }
    assert(cgnsFilePtr >= 0);
  }

  void DatabaseIO::closeDatabase() const
  {
    if (cgnsFilePtr != -1) {
      cg_close(cgnsFilePtr);
    }
    cgnsFilePtr = -1;
  }

  int64_t DatabaseIO::node_global_to_local(int64_t global, bool must_exist) const
  {
    return global;
  }

  int64_t DatabaseIO::element_global_to_local(int64_t global) const
  {
    return global;
  }

  void DatabaseIO::read_meta_data()
  {
    openDatabase();
    
    // Determine the number of bases in the grid.
    // Currently only handle 1.
    cgsize_t n_bases = 0;
    cg_nbases(cgnsFilePtr, &n_bases);
    if (n_bases != 1) {
      std::ostringstream errmsg;
      errmsg << "CGNS: Too many bases; only support files with a single bases at this time";
      IOSS_ERROR(errmsg);
    }

    // ========================================================================
    // Get the number of families in the mesh...
    // Will treat these as sidesets if they are of the type "FamilyBC_t"
    cgsize_t base = 1;
    cgsize_t num_families = 0;
    cg_nfamilies(cgnsFilePtr, base, &num_families);
    for (cgsize_t family=1; family <= num_families; family++) {
      char name[33];
      cgsize_t num_bc = 0;
      cgsize_t num_geo = 0;
      cg_family_read(cgnsFilePtr, base, family, name, &num_bc, &num_geo);
      std::cout << "Family " << family << " named " << name
		<< " has " << num_bc << " BC, and "
		<< num_geo << " geometry references\n";
      if (num_bc > 0) {
	// Create a sideset...
	std::string ss_name(name);
	Ioss::SideSet *ss = new Ioss::SideSet(this, ss_name);
	get_region()->add(ss);
      }
    }

    // ========================================================================
    // Get the number of zones (element blocks) in the mesh...
    int num_zones = 0;
    cg_nzones(cgnsFilePtr, base, &num_zones);
    m_blockLocalNodeMap.resize(num_zones+1);  // Let's use 1-based zones...
    m_zoneOffset.resize(num_zones+1);  // Let's use 1-based zones...
    
    // ========================================================================
    cgsize_t num_node = 0;
    cgsize_t num_elem = 0;
    for (cgsize_t zone=1; zone <= num_zones; zone++) {
      CG_ZoneType_t zone_type;
      cg_zone_type(cgnsFilePtr, base, zone, &zone_type);

      // See if all zones are "Unstructured" which is all we currently support...
      if (zone_type != CG_Unstructured) {
        std::ostringstream errmsg;
        errmsg << "ERROR: CGNS: Zone " << zone
	       << " is not of type Unstructured which is the only type currently supported";
	IOSS_ERROR(errmsg);
      }
      else {
	cgsize_t size[3];
	char zone_name[33];
	cg_zone_read(cgnsFilePtr, base, zone, zone_name, size);
	m_zoneNameMap[zone_name] = zone;
	
	cgsize_t total_block_nodes = size[0];
	m_blockLocalNodeMap[zone].resize(total_block_nodes, -1);
	
	// Determine number of "shared" nodes (shared with other zones)
	int nconn = 0;
	cg_nconns(cgnsFilePtr, base, zone, &nconn);
	cgsize_t num_shared = 0;
	for (int i=0; i < nconn; i++) {
	  char connectname[33];
	  CG_GridLocation_t location;
	  CG_GridConnectivityType_t connect_type;
	  CG_PointSetType_t ptset_type;
	  cgsize_t npnts = 0;
	  char donorname[33];
	  CG_ZoneType_t donor_zonetype;
	  CG_PointSetType_t donor_ptset_type;
	  CG_DataType_t donor_datatype;
	  cgsize_t ndata_donor;
	  
	  cg_conn_info(cgnsFilePtr, base, zone, i+1, connectname,
		       &location, &connect_type,
		       &ptset_type, &npnts, donorname,
		       &donor_zonetype, &donor_ptset_type,
		       &donor_datatype, &ndata_donor);

	  if (connect_type != CG_Abutting1to1 ||
	      ptset_type != CG_PointList ||
	      donor_ptset_type != CG_PointListDonor) {
	    std::ostringstream errmsg;
	    errmsg << "ERROR: CGNS: Zone " << zone
		   << " adjacency data is not correct type. Require Abutting1to1 and PointList."
		   << connect_type << "\t" << ptset_type << "\t" << donor_ptset_type;
	    IOSS_ERROR(errmsg);
	  }

	  // Verify data consistency...
	  if (npnts != ndata_donor) {
	    std::ostringstream errmsg;
	    errmsg << "ERROR: CGNS: Zone " << zone
		   << " point count (" << npnts << ") does not match donor point count (" << ndata_donor << ").";
	    IOSS_ERROR(errmsg);
	  }

	  // Get number of nodes shared with other "previous" zones...
	  // A "previous" zone will have a lower zone number this this zone...
	  auto donor_iter = m_zoneNameMap.find(donorname);
	  if (donor_iter != m_zoneNameMap.end() && (*donor_iter).second < zone) {
	    num_shared += npnts;
	    std::cout << "Zone " << zone << " shares " << npnts << " nodes with " << donorname << "\n";

	    std::vector<cgsize_t> points(npnts);
	    std::vector<cgsize_t> donors(npnts);

	    cg_conn_read(cgnsFilePtr, base, zone, i+1, TOPTR(points),
			 donor_datatype, TOPTR(donors));
	    
	    // Fill in entries in m_blockLocalNodeMap for the shared nodes...
	    auto &donor_map = m_blockLocalNodeMap[(*donor_iter).second];
	    auto &block_map = m_blockLocalNodeMap[zone];
	    for (int j = 0; j < npnts; j++) {
	      cgsize_t point = points[j];
	      cgsize_t donor = donors[j];
	      block_map[point-1] = donor_map[donor-1];
	    }
	  }
	}
	
	auto &block_map = m_blockLocalNodeMap[zone];
	size_t offset = num_node;
	for (cgsize_t i=0; i < total_block_nodes; i++) {
	  if (block_map[i] == -1) {
	    block_map[i] = offset++;
	  }
	}
	num_node = offset;

	cgsize_t total_elements = size[1];
	m_zoneOffset[zone] = num_elem;
	num_elem += total_elements;
	
	// NOTE: A Zone will have a single set of nodes, but can have
	//       multiple sections each with their own element type...
	//       Keep treating sections as element blocks until we
	//       have handled 'size[1]' number of elements; the remaining
	//       sections are then the boundary faces (?)
	cgsize_t num_sections = 0;
	cg_nsections(cgnsFilePtr, base, zone, &num_sections);
	
	// ========================================================================
	// Read the sections and create an element block for the ones that
	// define elements.  Some define boundary conditions...
	Ioss::ElementBlock *eblock = nullptr;
	
	
	for (cgsize_t is = 1; is <= num_sections; is++) {
	  char section_name[33];
	  CG_ElementType_t e_type;
	  cgsize_t el_start = 0;
	  cgsize_t el_end = 0;
	  cgsize_t num_bndry = 0;
	  cgsize_t parent_flag = 0;

	  // Get the type of elements in this section...
	  cg_section_read(cgnsFilePtr, base, zone, is,
			  section_name, &e_type, &el_start, &el_end, &num_bndry, &parent_flag);

	  cgsize_t num_entity = el_end - el_start + 1;

	  if (parent_flag == 0 && total_elements > 0) {
	    total_elements -= num_entity;
	    std::string element_topo = map_cgns_to_topology_type(e_type);
	    std::cout << "Added block " << zone_name
		      << ": CGNS topology = '" << cg_ElementTypeName(e_type)
		      << "', IOSS topology = '" << element_topo
		      << "' with " << num_entity << " elements\n";

	    eblock = new Ioss::ElementBlock(this, zone_name, element_topo, num_entity);
	    eblock->property_add(Ioss::Property("base", base));
	    eblock->property_add(Ioss::Property("zone", zone));
	    eblock->property_add(Ioss::Property("section", is));

	    assert(is == 1); // For now, assume each zone has only a single element block.
	    bool added = get_region()->add(eblock);
	    if(!added) {
	      delete eblock;
	      eblock = nullptr;
	    }
	  }
	  else {
	    // This is a boundary-condition -- sideset (?)
	    // See if there is an existing sideset with this name...
	    Ioss::SideSet *sset = get_region()->get_sideset(section_name);

#if 0
	    // NOTE: This if block is assuming that all BC (sidesets) are listed
	    //       as Family nodes.  While iterating zones, only create a BC
	    //       if this already exists as a sideset from the family nodes.
	    if (sset == nullptr) {
	      sset = new Ioss::SideSet(this, section_name);
	      bool added = get_region()->add(sset);
	      if(!added) {
		std::cout << "ERROR: Could not add sideset " << section_name << "\n";
		delete sset;
		sset = nullptr;
	      }
	    }
#endif
	    if (sset != nullptr) {
	      std::string block_name(zone_name);
	      block_name += "/";
	      block_name += section_name;
	      std::string face_topo = map_cgns_to_topology_type(e_type);
	      std::cout << "Added sideset " << block_name << " of topo " << face_topo
			<< " with " << num_entity << " faces\n";
	      
	      std::string parent_topo = eblock == nullptr ? "unknown" : eblock->topology()->name();
	      Ioss::SideBlock *sblk = new Ioss::SideBlock(this, block_name, face_topo, parent_topo,
							  num_entity);
	      sblk->property_add(Ioss::Property("base", base));
	      sblk->property_add(Ioss::Property("zone", zone));
	      sblk->property_add(Ioss::Property("section", is));
	      if (eblock != nullptr) {
		sblk->set_parent_element_block(eblock);
	      }
	      sset->add(sblk);
	    }
	  }
	}
      }
    }
    Ioss::NodeBlock *nblock = new Ioss::NodeBlock(this, "nodeblock_1", num_node, 3);
    nblock->property_add(Ioss::Property("base", base));
    get_region()->add(nblock);
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

  int64_t DatabaseIO::get_field_internal(const Ioss::Region* /* reg */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock* nb, const Ioss::Field& field,
					 void *data, size_t data_size) const
  {
    Ioss::Field::RoleType role = field.get_role();
    cgsize_t base = nb->get_property("base").get_int();
    cgsize_t num_to_get = field.verify(data_size);
    cgsize_t first = 1;

    char basename[33];

    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "mesh_model_coordinates_x") {
	double *rdata = static_cast<double*>(data);

	for (int zone=1; zone < (int)m_blockLocalNodeMap.size(); zone++) {
	  auto &block_map = m_blockLocalNodeMap[zone];     
	  cgsize_t num_coord = block_map.size();
	  std::vector<double> coord(num_coord);
	  int ierr = cg_coord_read(cgnsFilePtr, base, zone, "CoordinateX", CG_RealDouble,
				   &first, &num_coord, TOPTR(coord));
	  if (ierr < 0) {
	    cgns_error(cgnsFilePtr, __LINE__, myProcessor);
	  }

	  // Map to global coordinate position...
	  for (cgsize_t i=0; i < num_coord; i++) {
	    rdata[block_map[i]] = coord[i];
	  }
	}
      }
      
      else if (field.get_name() == "mesh_model_coordinates_y") {
	double *rdata = static_cast<double*>(data);

	for (int zone=1; zone < (int)m_blockLocalNodeMap.size(); zone++) {
	  auto &block_map = m_blockLocalNodeMap[zone];     
	  cgsize_t num_coord = block_map.size();
	  std::vector<double> coord(num_coord);
	  int ierr = cg_coord_read(cgnsFilePtr, base, zone, "CoordinateY", CG_RealDouble,
				   &first, &num_coord, TOPTR(coord));
	  if (ierr < 0) {
	    cgns_error(cgnsFilePtr, __LINE__, myProcessor);
}

	  // Map to global coordinate position...
	  for (cgsize_t i=0; i < num_coord; i++) {
	    rdata[block_map[i]] = coord[i];
	  }
	}
      }
      
      else if (field.get_name() == "mesh_model_coordinates_z") {
	double *rdata = static_cast<double*>(data);

	for (int zone=1; zone < (int)m_blockLocalNodeMap.size(); zone++) {
	  auto &block_map = m_blockLocalNodeMap[zone];     
	  cgsize_t num_coord = block_map.size();
	  std::vector<double> coord(num_coord);
	  int ierr = cg_coord_read(cgnsFilePtr, base, zone, "CoordinateZ", CG_RealDouble,
				   &first, &num_coord, TOPTR(coord));
	  if (ierr < 0) {
	    cgns_error(cgnsFilePtr, __LINE__, myProcessor);
}

	  // Map to global coordinate position...
	  for (cgsize_t i=0; i < num_coord; i++) {
	    rdata[block_map[i]] = coord[i];
	  }
	}
      }

      else if (field.get_name() == "mesh_model_coordinates") {
	cgsize_t cell_dimension = 0;
	cgsize_t phys_dimension = 0;
	cg_base_read(cgnsFilePtr, base, basename, &cell_dimension, &phys_dimension);

	double *rdata = static_cast<double*>(data);

	// Data required by upper classes store x0, y0, z0, ... xn,
	// yn, zn. Data stored in exodusII file is x0, ..., xn, y0,
	// ..., yn, z0, ..., zn so we have to allocate some scratch
	// memory to read in the data and then map into supplied
	// 'data'
	for (int zone=1; zone < (int)m_blockLocalNodeMap.size(); zone++) {
	  auto &block_map = m_blockLocalNodeMap[zone];     
	  cgsize_t num_coord = block_map.size();
	  std::vector<double> coord(num_coord);

	  int ierr = cg_coord_read(cgnsFilePtr, base, zone, "CoordinateX", CG_RealDouble,
				   &first, &num_coord, TOPTR(coord));
	  if (ierr < 0) {
	    cgns_error(cgnsFilePtr, __LINE__, myProcessor);
}

	  // Map to global coordinate position...
	  for (cgsize_t i=0; i < num_coord; i++) {
	    rdata[phys_dimension*block_map[i]+0] = coord[i];
	  }

	  ierr = cg_coord_read(cgnsFilePtr, base, zone, "CoordinateY", CG_RealDouble,
			       &first, &num_coord, TOPTR(coord));
	  if (ierr < 0) {
	    cgns_error(cgnsFilePtr, __LINE__, myProcessor);
}

	  // Map to global coordinate position...
	  for (cgsize_t i=0; i < num_coord; i++) {
	    rdata[phys_dimension*block_map[i]+1] = coord[i];
	  }

	  ierr = cg_coord_read(cgnsFilePtr, base, zone, "CoordinateZ", CG_RealDouble,
			       &first, &num_coord, TOPTR(coord));
	  if (ierr < 0) {
	    cgns_error(cgnsFilePtr, __LINE__, myProcessor);
}

	  // Map to global coordinate position...
	  for (cgsize_t i=0; i < num_coord; i++) {
	    rdata[phys_dimension*block_map[i]+2] = coord[i];
	  }
	}
      }
      else if (field.get_name() == "ids") {
	// Map the local ids in this node block
	// (1...node_count) to global node ids.
	if (field.get_type() == Ioss::Field::INT64) {
	  int64_t *idata = static_cast<int64_t*>(data);
	  std::iota(idata, idata+num_to_get, 1);
	}
	else {
	  assert(field.get_type() == Ioss::Field::INT32);
	  int *idata = static_cast<int*>(data);
	  std::iota(idata, idata+num_to_get, 1);
	}
      }
      else {
	num_to_get = Ioss::Utils::field_warning(nb, field, "input");
      }
      return num_to_get;
    }
    return -1;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeBlock* /* nb */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceBlock* /* nb */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock* eb,
					 const Ioss::Field& field,
					 void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {

      cgsize_t base = eb->get_property("base").get_int();
      cgsize_t zone = eb->get_property("zone").get_int();
      cgsize_t sect = eb->get_property("section").get_int();
      size_t my_element_count = eb->get_property("entity_count").get_int();
      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
	// Handle the MESH fields required for a CGNS file model.
	// (The 'genesis' portion)

	if (field.get_name() == "connectivity") {
	  // TODO: Need to map local to global...
	  cgsize_t *idata = (cgsize_t*)data;
	  int element_nodes = eb->topology()->number_nodes();
	  assert(field.raw_storage()->component_count() == element_nodes);

	  if (my_element_count > 0) {
	    int ierr = cg_elements_read(cgnsFilePtr, base, zone, sect,
					idata, nullptr);
	    if (ierr < 0) {
	      cgns_error(cgnsFilePtr, __LINE__, myProcessor);
}
	  }

	  // Now need to map block-local node connectivity to global nodes...
	  const auto &block_map = m_blockLocalNodeMap[zone];
	  for (size_t i=0; i < element_nodes * num_to_get; i++) {
	    idata[i] = block_map[idata[i]-1]+1;
	  }
	}
	else if (field.get_name() == "connectivity_raw") {
	  assert(field.raw_storage()->component_count() == eb->topology()->number_nodes());

	  if (my_element_count > 0) {
	    int ierr = cg_elements_read(cgnsFilePtr, base, zone, sect,
					(cgsize_t*)data, nullptr);
	    if (ierr < 0) {
	      cgns_error(cgnsFilePtr, __LINE__, myProcessor);
}
	  }
	}
	else if (field.get_name() == "ids") {
	  // TODO: This needs to change for parallel.
	  // Map the local ids in this element block
	  // (eb_offset+1...eb_offset+1+my_element_count) to global element ids.
	  size_t eb_offset_plus_one = eb->get_offset() + 1;
	  if (field.get_type() == Ioss::Field::INT64) {
	    int64_t *idata = static_cast<int64_t*>(data);
	    std::iota(idata, idata+my_element_count, eb_offset_plus_one);
	  } else {
	    assert(field.get_type() == Ioss::Field::INT32);
	    int *idata = static_cast<int*>(data);
	    std::iota(idata, idata+my_element_count, eb_offset_plus_one);
	  }
	}
	else if (field.get_name() == "implicit_ids") {
	  // TODO: This needs to change for parallel.
	  // If not parallel, then this is just 
	  // (eb_offset+1...eb_offset+1+my_element_count).
	  size_t eb_offset_plus_one = eb->get_offset() + 1;
	  if (field.get_type() == Ioss::Field::INT64) {
	    int64_t *idata = static_cast<int64_t*>(data);
	    std::iota(idata, idata+my_element_count, eb_offset_plus_one);
	  } else {
	    assert(field.get_type() == Ioss::Field::INT32);
	    int *idata = static_cast<int*>(data);
	    std::iota(idata, idata+my_element_count, eb_offset_plus_one);
	  }
	}
	else {
	  num_to_get = Ioss::Utils::field_warning(eb, field, "input");
	}
      }
      else {
	num_to_get = Ioss::Utils::field_warning(eb, field, "unknown");
      }
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet* /* ns */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeSet* /* ns */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceSet* /* ns */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementSet* /* ns */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock* sb, const Ioss::Field& field,
					 void *data , size_t data_size) const
  {
    cgsize_t base = sb->get_property("base").get_int();
    cgsize_t zone = sb->get_property("zone").get_int();
    cgsize_t sect = sb->get_property("section").get_int();
    
    ssize_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {
      int64_t entity_count = sb->get_property("entity_count").get_int();
      if (num_to_get != entity_count) {
	std::ostringstream errmsg;
	errmsg << "ERROR: Partial field input not yet implemented for side blocks";
	IOSS_ERROR(errmsg);
      }
    }

    Ioss::Field::RoleType role = field.get_role();
    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "element_side_raw" || field.get_name() == "element_side") {

	// TODO? Possibly rewrite using cgi_read_int_data so can skip reading element connectivity
	int nodes_per_face = sb->topology()->number_nodes();
	std::vector<cgsize_t> elements(nodes_per_face*num_to_get); // Not needed, but can't skip

	// We get:
	// *  num_to_get parent elements,
	// *  num_to_get zeros (other parent element for face, but on boundary so 0)
	// *  num_to_get face_on_element
	// *  num_to_get zeros (face on other parent element)
	std::vector<cgsize_t> parent(4 * num_to_get);

	int ierr = cg_elements_read(cgnsFilePtr, base, zone, sect,
				    TOPTR(elements), TOPTR(parent));
	if (ierr < 0) {
	  cgns_error(cgnsFilePtr, __LINE__, myProcessor);
}

	size_t offset = m_zoneOffset[zone];
	if (field.get_type() == Ioss::Field::INT32) {
	  int *idata = (int*)data;
	  size_t j = 0;
	  for (ssize_t i=0; i < num_to_get; i++) {
	    idata[j++] = parent[num_to_get*0 + i]+offset;  // Element
	    idata[j++] = parent[num_to_get*2 + i];
	    assert(parent[num_to_get*1+i] == 0);
	    assert(parent[num_to_get*3+i] == 0);
	  }
	  // Adjust face numbers to IOSS convention instead of CGNS convention...
	  map_cgns_face_to_ioss(sb->parent_element_topology(), num_to_get, idata);
	}
	else {
	  int64_t *idata = (int64_t*)data;
	  size_t j = 0;
	  for (ssize_t i=0; i < num_to_get; i++) {
	    idata[j++] = parent[num_to_get*0 + i]+offset; // Element
	    idata[j++] = parent[num_to_get*2 + i];
	    assert(parent[num_to_get*1+i] == 0);
	    assert(parent[num_to_get*3+i] == 0);
	  }
	  // Adjust face numbers to IOSS convention instead of CGNS convention...
	  map_cgns_face_to_ioss(sb->parent_element_topology(), num_to_get, idata);
	}


      }
    }
    return -1;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::SideSet* /* fs */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::CommSet* /* cs */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::Region* region, const Ioss::Field& field,
					 void *data, size_t data_size) const
  {
    return -1;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock* /* eb */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock* /* nb */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock* /* nb */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock* /* nb */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet* /* ns */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeSet* /* ns */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::FaceSet* /* ns */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::ElementSet* /* ns */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock* /* fb */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::SideSet* /* fs */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::CommSet* /* cs */, const Ioss::Field& /* field */,
					 void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  unsigned DatabaseIO::entity_field_support() const
  {
    return Ioss::REGION;
  }
}

