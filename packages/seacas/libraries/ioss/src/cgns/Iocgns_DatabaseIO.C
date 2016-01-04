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
      errmsg << "ERROR: CGNS Currently only supports reading, not writing.\n";
      IOSS_ERROR(errmsg);
    }
    openDatabase();
  }

  DatabaseIO::~DatabaseIO()
  { }

  void DatabaseIO::openDatabase() const
  {
    if (cgnsFilePtr < 0) {
      if (is_input()) {
	std::string decoded_filename = util().decode_filename(get_filename(),
							      isParallel);
        int ierr = cg_open(decoded_filename.c_str(), CG_MODE_READ, &cgnsFilePtr);
	if (ierr != CG_OK) {
	  // NOTE: Code will not continue past this call...
	  std::ostringstream errmsg;
	  errmsg << "ERROR: Problem opening file '" << decoded_filename << "' for read access.";
	  IOSS_ERROR(errmsg);
	}
      }
    }

    assert(cgnsFilePtr >= 0);
  }

  void DatabaseIO::closeDatabase() const
  {
    if (cgnsFilePtr != -1)
      cg_close(cgnsFilePtr);
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
    if (cg_nbases(cgnsFilePtr, &n_bases) != CG_OK) {
    }
    if (n_bases != 1) {
      std::ostringstream errmsg;
      errmsg << "CGNS: Too many bases; only support files with a single bases at this time";
      IOSS_ERROR(errmsg);
    }

    // Get the cell/element and physical dimension...
    cgsize_t base = 1;
    char basename[33];
    cgsize_t cell_dimension = 0;
    cgsize_t phys_dimension = 0;
    cg_base_read(cgnsFilePtr, base, basename, &cell_dimension, &phys_dimension);
    
    // ========================================================================
    // Get the number of families in the mesh...
    // Will treat these as sidesets if they are of the type "FamilyBC_t"
    cgsize_t num_families = 0;
    cg_nfamilies(cgnsFilePtr, base, &num_families);
    for (cgsize_t family=1; family <= num_families; family++) {
      char name[33];
      cgsize_t num_bc = 0;
      cgsize_t num_geo = 0;
      cg_family_read(cgnsFilePtr, base, family, name, &num_bc, &num_geo);
      std::cerr << "Family " << family << " named " << name
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
    cgsize_t num_zones = 0;
    cg_nzones(cgnsFilePtr, base, &num_zones);

    
    // ========================================================================
    size_t num_node = 0;
    size_t num_elem = 0;
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

	{
	  std::string block_name(zone_name);
	  block_name += "_nodes";
	  Ioss::NodeBlock *block = new Ioss::NodeBlock(this, block_name, size[0], phys_dimension);
	  block->property_add(Ioss::Property("zone", zone));
	  block->property_add(Ioss::Property("base", base));
	  bool added = get_region()->add(block);
	  if(!added) {
	    delete block;
	  }
	}
	
	cgsize_t total_elements = size[1];
	
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
	Ioss::ElementBlock *eblock = NULL;
	
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
	    std::cerr << "Added block " << zone_name << " of topo " << element_topo
		      << " with " << num_entity << " elements\n";

	    eblock = new Ioss::ElementBlock(this, zone_name, element_topo, num_entity);
	    eblock->property_add(Ioss::Property("base", base));
	    eblock->property_add(Ioss::Property("zone", zone));
	    eblock->property_add(Ioss::Property("section", is));
	    bool added = get_region()->add(eblock);
	    if(!added) {
	      delete eblock;
	      eblock = NULL;
	    }
	  }
	  else {
	    // This is a boundary-condition -- sideset (?)
	    // See if there is an existing sideset with this name...
	    Ioss::SideSet *sset = get_region()->get_sideset(section_name);
	    if (sset == NULL) {
	      Ioss::SideSet *sset = new Ioss::SideSet(this, section_name);
	      bool added = get_region()->add(sset);
	      if(!added) {
		delete sset;
		sset = NULL;
	      }
	    }

	    if (sset != NULL) {
	      std::string block_name(zone_name);
	      block_name += "/";
	      block_name += section_name;
	      std::string face_topo = map_cgns_to_topology_type(e_type);
	      std::cerr << "Added sideset " << block_name << " of topo " << face_topo
			<< " with " << num_entity << " faces\n";
	      
	      std::string parent_topo = eblock == NULL ? "unknown" : eblock->topology()->name();
	      Ioss::SideBlock *sblk = new Ioss::SideBlock(this, block_name, face_topo, parent_topo,
							  num_entity);
	      sblk->property_add(Ioss::Property("base", base));
	      sblk->property_add(Ioss::Property("zone", zone));
	      sblk->property_add(Ioss::Property("section", is));
	      if (eblock != NULL) {
		sblk->set_parent_element_block(eblock);
	      }
	      sset->add(sblk);
	    }
	  }
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
    cgsize_t zone = nb->get_property("zone").get_int();

#if 0
    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "mesh_model_coordinates_x") {
	double *rdata = static_cast<double*>(data);


	int ierr = ex_get_coord(get_file_pointer(), rdata, NULL, NULL);
	if (ierr < 0)
	  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
      
      else if (field.get_name() == "mesh_model_coordinates_y") {
	double *rdata = static_cast<double*>(data);
	int ierr = ex_get_coord(get_file_pointer(), NULL, rdata, NULL);
	if (ierr < 0)
	  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }
      
      else if (field.get_name() == "mesh_model_coordinates_z") {
	double *rdata = static_cast<double*>(data);
	int ierr = ex_get_coord(get_file_pointer(), NULL, NULL, rdata);
	if (ierr < 0)
	  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);
      }

      else if (field.get_name() == "mesh_model_coordinates") {
	// Data required by upper classes store x0, y0, z0, ... xn,
	// yn, zn. Data stored in exodusII file is x0, ..., xn, y0,
	// ..., yn, z0, ..., zn so we have to allocate some scratch
	// memory to read in the data and then map into supplied
	// 'data'
	std::vector<double> x(num_to_get);
	std::vector<double> y;
	if (spatialDimension > 1)
	  y.resize(num_to_get);
	std::vector<double> z;
	if (spatialDimension == 3)
	  z.resize(num_to_get);

	// Cast 'data' to correct size -- double
	double *rdata = static_cast<double*>(data);

	int ierr = ex_get_coord(get_file_pointer(), TOPTR(x), TOPTR(y), TOPTR(z));
	if (ierr < 0)
	  Ioex::exodus_error(get_file_pointer(), __LINE__, myProcessor);

	size_t index = 0;
	for (size_t i=0; i < num_to_get; i++) {
	  rdata[index++] = x[i];
	  if (spatialDimension > 1)
	    rdata[index++] = y[i];
	  if (spatialDimension == 3)
	    rdata[index++] = z[i];
	}
      }

      else if (field.get_name() == "ids") {
	// Map the local ids in this node block
	// (1...node_count) to global node ids.
	get_map(EX_NODE_BLOCK).map_implicit_data(data, field, num_to_get, 0);
      }
#endif
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
    int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock* /* eb */, const Ioss::Field& /* field */,
					   void */* data */, size_t /* data_size */) const
    {
      return -1;
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
    int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock* /* eb */, const Ioss::Field& /* field */,
					   void */* data */, size_t /* data_size */) const
    {
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

