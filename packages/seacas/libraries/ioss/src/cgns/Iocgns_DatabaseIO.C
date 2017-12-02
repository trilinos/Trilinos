// CGNS Assumptions:
// * All boundary conditions are listed as Family nodes at the "top" level.
// * Single element block per zone.
// * Single Base.
// * ZoneGridConnectivity is 1to1 with point lists for unstructured

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

#include <Ioss_CodeTypes.h>
#include <Ioss_Utils.h>
#include <cassert>
#include <cgns/Iocgns_DatabaseIO.h>
#include <cgns/Iocgns_Utils.h>
#include <cgnslib.h>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <sys/select.h>
#include <vector>

#if !defined(CGNSLIB_H)
#error "Could not include cgnslib.h"
#endif

#include "Ioss_CommSet.h"
#include "Ioss_DBUsage.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_EdgeBlock.h"
#include "Ioss_EdgeSet.h"
#include "Ioss_ElementBlock.h"
#include "Ioss_ElementSet.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_EntityType.h"
#include "Ioss_FaceBlock.h"
#include "Ioss_FaceSet.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_NodeSet.h"
#include "Ioss_SideBlock.h"
#include "Ioss_SideSet.h"
#include "Ioss_StructuredBlock.h"

#include "Ioss_Field.h"
#include "Ioss_IOFactory.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_Property.h"
#include "Ioss_Region.h"
#include "Ioss_State.h"
#include "Ioss_Utils.h"
#include "Ioss_VariableType.h"

extern char hdf5_access[64];

namespace Iocgns {

  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string &filename,
                         Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
                         const Ioss::PropertyManager &props)
      : Ioss::DatabaseIO(region, filename, db_usage, communicator, props), cgnsFilePtr(-1),
        nodeCount(0), elementCount(0)
  {
    dbState = Ioss::STATE_UNKNOWN;

#if IOSS_DEBUG_OUTPUT
    std::cout << "CGNS DatabaseIO using " << CG_SIZEOF_SIZE << "-bit integers.\n";
#endif
    openDatabase__();
  }

  DatabaseIO::~DatabaseIO()
  {
    for (auto &gtb : m_globalToBlockLocalNodeMap) {
      delete gtb.second;
    }
    cg_close(cgnsFilePtr);
  }

  void DatabaseIO::openDatabase__() const
  {
    if (cgnsFilePtr < 0) {
      if ((is_input() && properties.exists("MEMORY_READ")) ||
          (!is_input() && properties.exists("MEMORY_WRITE"))) {
        strcpy(hdf5_access, "PARALLEL");
      }
      int mode = is_input() ? CG_MODE_READ : CG_MODE_WRITE;
      int ierr = cg_open(get_filename().c_str(), mode, &cgnsFilePtr);
      if (ierr != CG_OK) {
        // NOTE: Code will not continue past this call...
        std::ostringstream errmsg;
        errmsg << "ERROR: Problem opening file '" << get_filename() << "' for "
               << (is_input() ? "read" : "write") << " access. "
               << "CGNS Error: '" << cg_get_error() << "'";
        IOSS_ERROR(errmsg);
      }
      if ((is_input() && properties.exists("MEMORY_READ")) ||
          (!is_input() && properties.exists("MEMORY_WRITE"))) {
        strcpy(hdf5_access, "NATIVE");
      }

      if (properties.exists("INTEGER_SIZE_API")) {
        int isize = properties.get("INTEGER_SIZE_API").get_int();
        if (isize == 8) {
          set_int_byte_size_api(Ioss::USE_INT64_API);
        }
        if (isize == 4) {
          set_int_byte_size_api(Ioss::USE_INT32_API);
        }
      }
      else if (CG_SIZEOF_SIZE == 64) {
        set_int_byte_size_api(Ioss::USE_INT64_API);
      }

#if 0
      // This isn't currently working since CGNS currently has chunking
      // disabled for HDF5 files and compression requires chunking.
      if (!is_input()) {
	if (properties.exists("COMPRESSION_LEVEL")) {
	  int comp = properties.get("COMPRESSION_LEVEL").get_int();
	  cg_configure(CG_CONFIG_HDF5_COMPRESS, (void*)comp);
	}
      }
#endif
    }
    assert(cgnsFilePtr >= 0);
  }

  void DatabaseIO::closeDatabase__() const
  {
    if (cgnsFilePtr != -1) {
      CGCHECK(cg_close(cgnsFilePtr));
    }
    cgnsFilePtr = -1;
  }

  void DatabaseIO::finalize_database()
  {
    if (is_input()) {
      return;
    }

    if (m_timesteps.empty()) {
      return;
    }

    Utils::finalize_database(cgnsFilePtr, m_timesteps, get_region(), myProcessor);
  }

  int64_t DatabaseIO::node_global_to_local__(int64_t global, bool /*must_exist*/) const
  {
    return global;
  }

  int64_t DatabaseIO::element_global_to_local__(int64_t global) const { return global; }

  void DatabaseIO::create_structured_block(int base, int zone, size_t &num_node)
  {
    cgsize_t size[9];
    char     zone_name[33];
    CGCHECK(cg_zone_read(cgnsFilePtr, base, zone, zone_name, size));
    m_zoneNameMap[zone_name] = zone;

    assert(size[0] - 1 == size[3]);
    assert(size[1] - 1 == size[4]);
    assert(size[2] - 1 == size[5]);

    assert(size[6] == 0);
    assert(size[7] == 0);
    assert(size[8] == 0);

    int index_dim = 0;
    CGCHECK(cg_index_dim(cgnsFilePtr, base, zone, &index_dim));
    // An Ioss::StructuredBlock corresponds to a CG_Structured zone...
    Ioss::StructuredBlock *block =
        new Ioss::StructuredBlock(this, zone_name, index_dim, size[3], size[4], size[5]);

    block->property_add(Ioss::Property("base", base));
    block->property_add(Ioss::Property("zone", zone));
    block->property_add(Ioss::Property("id", zone));
    block->property_add(Ioss::Property("guid", zone));
    get_region()->add(block);

    num_node += block->get_property("node_count").get_int();

    // Handle zone-grid-connectivity...
    int nconn = 0;
    CGCHECK(cg_n1to1(cgnsFilePtr, base, zone, &nconn));
    for (int i = 0; i < nconn; i++) {
      char                    connectname[33];
      char                    donorname[33];
      std::array<cgsize_t, 6> range;
      std::array<cgsize_t, 6> donor_range;
      Ioss::IJK_t             transform;

      CGCHECK(cg_1to1_read(cgnsFilePtr, base, zone, i + 1, connectname, donorname, range.data(),
                           donor_range.data(), transform.data()));

      // Get number of nodes shared with other "previous" zones...
      // A "previous" zone will have a lower zone number this this zone...
      int  donor_zone = -1;
      auto donor_iter = m_zoneNameMap.find(donorname);
      if (donor_iter != m_zoneNameMap.end()) {
        donor_zone = (*donor_iter).second;
      }
      Ioss::IJK_t range_beg{{(int)range[0], (int)range[1], (int)range[2]}};
      Ioss::IJK_t range_end{{(int)range[3], (int)range[4], (int)range[5]}};
      Ioss::IJK_t donor_beg{{(int)donor_range[0], (int)donor_range[1], (int)donor_range[2]}};
      Ioss::IJK_t donor_end{{(int)donor_range[3], (int)donor_range[4], (int)donor_range[5]}};

      bool owns_nodes = zone < donor_zone || donor_zone == -1;
      block->m_zoneConnectivity.emplace_back(connectname, zone, donorname, donor_zone, transform,
                                             range_beg, range_end, donor_beg, donor_end,
                                             owns_nodes);
      block->m_zoneConnectivity.back().m_donorProcessor = 0;
      block->m_zoneConnectivity.back().m_ownerProcessor = 0;
    }

    // Handle boundary conditions...
    Utils::add_structured_boundary_conditions(cgnsFilePtr, block);
  }

  size_t DatabaseIO::finalize_structured_blocks()
  {
    const auto &blocks = get_region()->get_structured_blocks();

    // If there are any Structured blocks, need to iterate them and their 1-to-1 connections
    // and update the donor_zone id for zones that had not yet been processed at the time of
    // definition...
    for (auto &block : blocks) {
      for (auto &conn : block->m_zoneConnectivity) {
        if (conn.m_donorZone < 0) {
          auto donor_iter = m_zoneNameMap.find(conn.m_donorName);
          assert(donor_iter != m_zoneNameMap.end());
          conn.m_donorZone = (*donor_iter).second;
        }
        conn.m_donorGUID = conn.m_donorZone;
        conn.m_ownerGUID = conn.m_ownerZone;
      }
    }

    size_t num_nodes = Utils::resolve_nodes(*get_region(), myProcessor, false);
    return num_nodes;
  }

  void DatabaseIO::create_unstructured_block(int base, int zone, size_t &num_node)
  {
    cgsize_t size[9];
    char     zone_name[33];
    CGCHECK(cg_zone_read(cgnsFilePtr, base, zone, zone_name, size));
    m_zoneNameMap[zone_name] = zone;

    size_t total_block_nodes = size[0];
    m_blockLocalNodeMap[zone].resize(total_block_nodes, -1);

    // Determine number of "shared" nodes (shared with other zones)
    if (zone > 1) { // Donor zone is always lower numbered, so zone 1 has no donor zone.
      int nconn = 0;
      CGCHECK(cg_nconns(cgnsFilePtr, base, zone, &nconn));
      cgsize_t num_shared = 0;
      for (int i = 0; i < nconn; i++) {
        char                      connectname[33];
        CG_GridLocation_t         location;
        CG_GridConnectivityType_t connect_type;
        CG_PointSetType_t         ptset_type;
        cgsize_t                  npnts = 0;
        char                      donorname[33];
        CG_ZoneType_t             donor_zonetype;
        CG_PointSetType_t         donor_ptset_type;
        CG_DataType_t             donor_datatype;
        cgsize_t                  ndata_donor;

        CGCHECK(cg_conn_info(cgnsFilePtr, base, zone, i + 1, connectname, &location, &connect_type,
                             &ptset_type, &npnts, donorname, &donor_zonetype, &donor_ptset_type,
                             &donor_datatype, &ndata_donor));

        if (connect_type != CG_Abutting1to1 || ptset_type != CG_PointList ||
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
          errmsg << "ERROR: CGNS: Zone " << zone << " point count (" << npnts
                 << ") does not match donor point count (" << ndata_donor << ").";
          IOSS_ERROR(errmsg);
        }

        // Get number of nodes shared with other "previous" zones...
        // A "previous" zone will have a lower zone number this this zone...
        auto donor_iter = m_zoneNameMap.find(donorname);
        if (donor_iter != m_zoneNameMap.end() && (*donor_iter).second < zone) {
          num_shared += npnts;
#if IOSS_DEBUG_OUTPUT
          std::cout << "Zone " << zone << " shares " << npnts << " nodes with " << donorname
                    << "\n";
#endif
          std::vector<cgsize_t> points(npnts);
          std::vector<cgsize_t> donors(npnts);

          CGCHECK(cg_conn_read(cgnsFilePtr, base, zone, i + 1, TOPTR(points), donor_datatype,
                               TOPTR(donors)));

          // Fill in entries in m_blockLocalNodeMap for the shared nodes...
          auto &donor_map = m_blockLocalNodeMap[(*donor_iter).second];
          auto &block_map = m_blockLocalNodeMap[zone];
          for (int j = 0; j < npnts; j++) {
            cgsize_t point       = points[j];
            cgsize_t donor       = donors[j];
            block_map[point - 1] = donor_map[donor - 1];
          }
        }
      }
    }

    auto & block_map = m_blockLocalNodeMap[zone];
    size_t offset    = num_node;
    for (size_t i = 0; i < total_block_nodes; i++) {
      if (block_map[i] == -1) {
        block_map[i] = offset++;
      }
    }
    num_node = offset;

    size_t num_elem = size[1];
    m_zoneOffset[zone]    = m_zoneOffset[zone - 1] + num_elem;

    // NOTE: A Zone will have a single set of nodes, but can have
    //       multiple sections each with their own element type...
    //       Keep treating sections as element blocks until we
    //       have handled 'size[1]' number of elements; the remaining
    //       sections are then the boundary faces (?)
    int num_sections = 0;
    CGCHECK(cg_nsections(cgnsFilePtr, base, zone, &num_sections));

    // ========================================================================
    // Read the sections and create an element block for the ones that
    // define elements.  Some define boundary conditions...
    Ioss::ElementBlock *eblock = nullptr;

    for (int is = 1; is <= num_sections; is++) {
      char             section_name[33];
      CG_ElementType_t e_type;
      cgsize_t         el_start    = 0;
      cgsize_t         el_end      = 0;
      int              num_bndry   = 0;
      int              parent_flag = 0;

      // Get the type of elements in this section...
      CGCHECK(cg_section_read(cgnsFilePtr, base, zone, is, section_name, &e_type, &el_start,
                              &el_end, &num_bndry, &parent_flag));

      cgsize_t num_entity = el_end - el_start + 1;

      if (parent_flag == 0 && num_elem > 0) {
        num_elem -= num_entity;
        std::string element_topo = Utils::map_cgns_to_topology_type(e_type);
#if IOSS_DEBUG_OUTPUT
        std::cout << "Added block " << zone_name << ": CGNS topology = '"
                  << cg_ElementTypeName(e_type) << "', IOSS topology = '" << element_topo
                  << "' with " << num_entity << " elements\n";
#endif
        eblock = new Ioss::ElementBlock(this, zone_name, element_topo, num_entity);
        eblock->property_add(Ioss::Property("base", base));
        eblock->property_add(Ioss::Property("zone", zone));
        eblock->property_add(Ioss::Property("id", zone));
        eblock->property_add(Ioss::Property("guid", zone));
        eblock->property_add(Ioss::Property("section", is));
        eblock->property_add(Ioss::Property("node_count", (int64_t)total_block_nodes));

        assert(is == 1); // For now, assume each zone has only a single element block.
        bool added = get_region()->add(eblock);
        if (!added) {
          delete eblock;
          eblock = nullptr;
        }
      }
      else {
        // This is a boundary-condition -- sideset (?)
        // See if there is an existing sideset with this name...
        Ioss::SideSet *sset = get_region()->get_sideset(section_name);

        if (sset != nullptr) {
          std::string block_name(zone_name);
          block_name += "/";
          block_name += section_name;
          std::string face_topo = Utils::map_cgns_to_topology_type(e_type);
#if IOSS_DEBUG_OUTPUT
          std::cout << "Added sideset " << block_name << " of topo " << face_topo << " with "
                    << num_entity << " faces\n";
#endif
          std::string parent_topo = eblock == nullptr ? "unknown" : eblock->topology()->name();
          auto sblk = new Ioss::SideBlock(this, block_name, face_topo, parent_topo, num_entity);
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

  void DatabaseIO::read_meta_data__()
  {
    openDatabase__();

    // Determine the number of bases in the grid.
    // Currently only handle 1.
    int n_bases = 0;
    CGCHECK(cg_nbases(cgnsFilePtr, &n_bases));
    if (n_bases != 1) {
      std::ostringstream errmsg;
      errmsg << "CGNS: Too many bases; only support files with a single bases at this time";
      IOSS_ERROR(errmsg);
    }

    get_step_times__();

    // ========================================================================
    // Get the number of families in the mesh...
    // Will treat these as sidesets if they are of the type "FamilyBC_t"
    Utils::add_sidesets(cgnsFilePtr, this);

    // ========================================================================
    // Get the number of zones (element blocks) in the mesh...
    int num_zones = 0;
    int base      = 1;
    CGCHECK(cg_nzones(cgnsFilePtr, base, &num_zones));
    m_blockLocalNodeMap.resize(num_zones + 1); // Let's use 1-based zones...
    m_zoneOffset.resize(num_zones + 1);        // Let's use 1-based zones...

    // ========================================================================
    size_t        num_node         = 0;
    CG_ZoneType_t common_zone_type = Utils::check_zone_type(cgnsFilePtr);

    for (int zone = 1; zone <= num_zones; zone++) {
      if (common_zone_type == CG_Structured) {
        create_structured_block(base, zone, num_node);
      }
      else if (common_zone_type == CG_Unstructured) {
        create_unstructured_block(base, zone, num_node);
      }
      else {
        // This should be handled already in check_zone_type...
        std::ostringstream errmsg;
        errmsg << "ERROR: CGNS: Zone " << zone
               << " is not of type Unstructured or Structured "
                  "which are the only types currently supported";
        IOSS_ERROR(errmsg);
      }
    }

    if (common_zone_type == CG_Structured) {
      num_node = finalize_structured_blocks();
    }

    char basename[33];
    int  cell_dimension = 0;
    int  phys_dimension = 0;
    CGCHECK(cg_base_read(cgnsFilePtr, base, basename, &cell_dimension, &phys_dimension));
#if IOSS_DEBUG_OUTPUT
    std::cout << "Physical dimension = " << phys_dimension << "\n";
#endif

    Ioss::NodeBlock *nblock = new Ioss::NodeBlock(this, "nodeblock_1", num_node, phys_dimension);
    nblock->property_add(Ioss::Property("base", base));
    get_region()->add(nblock);

    Utils::add_transient_variables(cgnsFilePtr, m_timesteps, get_region(), myProcessor);
  }

  void DatabaseIO::write_meta_data()
  {
    int num_zones = get_region()->get_property("element_block_count").get_int() +
                    get_region()->get_property("structured_block_count").get_int();
    m_bcOffset.resize(num_zones + 1);   // use 1-based zones...
    m_zoneOffset.resize(num_zones + 1); // use 1-based zones...

    Utils::common_write_meta_data(cgnsFilePtr, *get_region(), m_zoneOffset);
  }

  void DatabaseIO::get_step_times__()
  {
    Utils::get_step_times(cgnsFilePtr, m_timesteps, get_region(), timeScaleFactor, myProcessor);
  }

  void DatabaseIO::write_adjacency_data()
  {
    // Determine adjacency information between unstructured blocks.
    // Could save this information from the input mesh, but then
    // could not read an exodus mesh and write a cgns mesh.
    // However, in long run may still want to read/save input adjacency
    // data if doing cgns -> cgns...  For now, try generating information.

    // If block I is adjacent to block J, then they will share at
    // least 1 "side" (face 3D or edge 2D).
    // Currently, assuming they are adjacent if they share at least one node...

    size_t node_count = get_region()->get_property("node_count").get_int();

    const auto &blocks = get_region()->get_element_blocks();
    for (auto I = blocks.cbegin(); I != blocks.cend(); I++) {
      int base = (*I)->get_property("base").get_int();
      int zone = (*I)->get_property("zone").get_int();

      const auto &I_map = m_globalToBlockLocalNodeMap[zone];

      // Flag all nodes used by this block...
      std::vector<size_t> I_nodes(node_count);
      for (size_t i = 0; i < I_map->map().size() - 1; i++) {
        auto global     = I_map->map()[i + 1] - 1;
        I_nodes[global] = i + 1;
      }
      for (auto J = I + 1; J != blocks.end(); J++) {
        int                   dzone = (*J)->get_property("zone").get_int();
        const auto &          J_map = m_globalToBlockLocalNodeMap[dzone];
        std::vector<cgsize_t> point_list;
        std::vector<cgsize_t> point_list_donor;
        for (size_t i = 0; i < J_map->map().size() - 1; i++) {
          auto global = J_map->map()[i + 1] - 1;
          if (I_nodes[global] > 0) {
            // Have a match between nodes used by two different blocks,
            // They are adjacent...
            point_list.push_back(I_nodes[global]);
            point_list_donor.push_back(i + 1);
          }
        }

        // If point_list non_empty, then output this adjacency node...
        if (!point_list.empty()) {
          int         gc_idx = 0;
          std::string name   = (*I)->name();
          name += "_to_";
          name += (*J)->name();
          const auto &d1_name = (*J)->name();
          CGCHECK(cg_conn_write(cgnsFilePtr, base, zone, name.c_str(), CG_Vertex, CG_Abutting1to1,
                                CG_PointList, point_list.size(), TOPTR(point_list), d1_name.c_str(),
                                CG_Unstructured, CG_PointListDonor, CG_DataTypeNull,
                                point_list_donor.size(), TOPTR(point_list_donor), &gc_idx));

          name = (*J)->name();
          name += "_to_";
          name += (*I)->name();
          const auto &d2_name = (*I)->name();

          CGCHECK(cg_conn_write(cgnsFilePtr, base, dzone, name.c_str(), CG_Vertex, CG_Abutting1to1,
                                CG_PointList, point_list_donor.size(), TOPTR(point_list_donor),
                                d2_name.c_str(), CG_Unstructured, CG_PointListDonor,
                                CG_DataTypeNull, point_list.size(), TOPTR(point_list), &gc_idx));
        }
      }
    }
  }

  bool DatabaseIO::begin__(Ioss::State state)
  {
    dbState = state;
    return true;
  }

  bool DatabaseIO::end__(Ioss::State state)
  {
    // Transitioning out of state 'state'
    switch (state) {
    case Ioss::STATE_DEFINE_MODEL:
      if (!is_input() && open_create_behavior() != Ioss::DB_APPEND) {
        write_meta_data();
      }
      break;
    case Ioss::STATE_MODEL:
      if (!is_input() && open_create_behavior() != Ioss::DB_APPEND) {
        write_adjacency_data();
      }
      break;
    case Ioss::STATE_DEFINE_TRANSIENT:
      if (!is_input() && open_create_behavior() != Ioss::DB_APPEND) {
        write_results_meta_data();
      }
      break;
    default: // ignore everything else...
      break;
    }

    dbState = Ioss::STATE_UNKNOWN;
    return true;
  }

  bool DatabaseIO::begin_state__(Ioss::Region *region, int state, double time)
  {
    if (is_input()) {
      return true;
    }
    Utils::write_flow_solution_metadata(cgnsFilePtr, get_region(), state,
                                        &m_currentVertexSolutionIndex,
                                        &m_currentCellCenterSolutionIndex);

    return true;
  }

  bool DatabaseIO::end_state__(Ioss::Region * /* region */, int state, double time)
  {
    if (!is_input()) {
      m_timesteps.push_back(time);
      assert(m_timesteps.size() == (size_t)state);
    }
    return true;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::Region *reg, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(reg, field, "input");
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    Ioss::Field::RoleType role       = field.get_role();
    int                   base       = nb->get_property("base").get_int();
    size_t                num_to_get = field.verify(data_size);
    cgsize_t              first      = 1;

    char basename[33];

    // Create a lambda to eliminate lots of duplicate code in coordinate outputs...
    auto coord_lambda = [this, &data, &first, base](const char *ordinate) {
      double *rdata = static_cast<double *>(data);

      for (int zone = 1; zone < static_cast<int>(m_blockLocalNodeMap.size()); zone++) {
        auto &              block_map = m_blockLocalNodeMap[zone];
        cgsize_t            num_coord = block_map.size();
        std::vector<double> coord(num_coord);
        CGCHECK(cg_coord_read(cgnsFilePtr, base, zone, ordinate, CG_RealDouble, &first, &num_coord,
                              TOPTR(coord)));

        // Map to global coordinate position...
        for (cgsize_t i = 0; i < num_coord; i++) {
          rdata[block_map[i]] = coord[i];
        }
      }
    };
    // End of lambda...

    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "mesh_model_coordinates_x") {
        // Use the lambda...
        coord_lambda("CoordinateX");
      }

      else if (field.get_name() == "mesh_model_coordinates_y") {
        coord_lambda("CoordinateY");
      }

      else if (field.get_name() == "mesh_model_coordinates_z") {
        coord_lambda("CoordinateZ");
      }

      else if (field.get_name() == "mesh_model_coordinates") {
        int cell_dimension = 0;
        int phys_dimension = 0;
        CGCHECK(cg_base_read(cgnsFilePtr, base, basename, &cell_dimension, &phys_dimension));

        double *rdata = static_cast<double *>(data);

        // Data required by upper classes store x0, y0, z0, ... xn,
        // yn, zn. Data stored in exodusII file is x0, ..., xn, y0,
        // ..., yn, z0, ..., zn so we have to allocate some scratch
        // memory to read in the data and then map into supplied
        // 'data'
        for (int zone = 1; zone < static_cast<int>(m_blockLocalNodeMap.size()); zone++) {
          auto &              block_map = m_blockLocalNodeMap[zone];
          cgsize_t            num_coord = block_map.size();
          std::vector<double> coord(num_coord);

          // ========================================================================
          // Repetitive code for each coordinate direction; use a lambda to consolidate...
          auto blk_coord_lambda = [this, block_map, base, zone, &coord, first, num_coord,
                                   phys_dimension, &rdata](const char *ord_name, int ordinate) {
            CGCHECK(cg_coord_read(cgnsFilePtr, base, zone, ord_name, CG_RealDouble, &first,
                                  &num_coord, coord.data()));

            // Map to global coordinate position...
            for (cgsize_t i = 0; i < num_coord; i++) {
              rdata[phys_dimension * block_map[i] + ordinate] = coord[i];
            }
          };
          // End of lambda...
          // ========================================================================

          blk_coord_lambda("CoordinateX", 0);

          if (phys_dimension >= 2) {
            blk_coord_lambda("CoordinateY", 1);
          }

          if (phys_dimension >= 3) {
            blk_coord_lambda("CoordinateZ", 2);
          }
        }
      }
      else if (field.get_name() == "ids") {
        // Map the local ids in this node block
        // (1...node_count) to global node ids.
        if (field.get_type() == Ioss::Field::INT64) {
          int64_t *idata = static_cast<int64_t *>(data);
          std::iota(idata, idata + num_to_get, 1);
        }
        else {
          assert(field.get_type() == Ioss::Field::INT32);
          int *idata = static_cast<int *>(data);
          std::iota(idata, idata + num_to_get, 1);
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(nb, field, "input");
      }
    }
    else if (role == Ioss::Field::TRANSIENT) {
      // Locate the FlowSolution node corresponding to the correct state/step/time
      // TODO: do this at read_meta_data() and store...
      int step = get_region()->get_current_state();

      for (int zone = 1; zone < static_cast<int>(m_blockLocalNodeMap.size()); zone++) {
        int   solution_index = Utils::find_solution_index(cgnsFilePtr, base, zone, step, CG_Vertex);
        auto &block_map      = m_blockLocalNodeMap[zone];
        cgsize_t num_block_node = block_map.size();

        double *            rdata        = static_cast<double *>(data);
        cgsize_t            range_min[1] = {1};
        cgsize_t            range_max[1] = {num_block_node};
        auto                var_type     = field.transformed_storage();
        int                 comp_count   = var_type->component_count();
        std::vector<double> cgns_data(num_block_node);
        if (comp_count == 1) {
          CGCHECK(cg_field_read(cgnsFilePtr, base, zone, solution_index, field.get_name().c_str(),
                                CG_RealDouble, range_min, range_max, cgns_data.data()));

          // Map to global nodal field position...
          for (cgsize_t i = 0; i < num_block_node; i++) {
            rdata[block_map[i]] = cgns_data[i];
          }
        }
        else {
          char field_suffix_separator = get_field_separator();
          for (int i = 0; i < comp_count; i++) {
            std::string var_name =
                var_type->label_name(field.get_name(), i + 1, field_suffix_separator);

            CGCHECK(cg_field_read(cgnsFilePtr, base, zone, solution_index, var_name.c_str(),
                                  CG_RealDouble, range_min, range_max, cgns_data.data()));
            for (cgsize_t j = 0; j < num_block_node; j++) {
              auto global                    = block_map[j];
              rdata[comp_count * global + i] = cgns_data[j];
            }
          }
        }
      }
    }
    else {
      num_to_get = Ioss::Utils::field_warning(nb, field, "input");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(eb, field, "input");
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(fb, field, "input");
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {

      int                   base             = eb->get_property("base").get_int();
      int                   zone             = eb->get_property("zone").get_int();
      int                   sect             = eb->get_property("section").get_int();
      cgsize_t              my_element_count = eb->entity_count();
      Ioss::Field::RoleType role             = field.get_role();

      if (role == Ioss::Field::MESH) {
        // Handle the MESH fields required for a CGNS file model.
        // (The 'genesis' portion)

        if (field.get_name() == "connectivity" || field.get_name() == "connectivity_raw") {
          // TODO(gdsjaar): Need to map local to global...
          int element_nodes = eb->topology()->number_nodes();
          assert(field.raw_storage()->component_count() == element_nodes);

          if (my_element_count > 0) {
            int field_byte_size = (field.get_type() == Ioss::Field::INT32) ? 32 : 64;
            if (field_byte_size == CG_SIZEOF_SIZE) {
              cgsize_t *idata = reinterpret_cast<cgsize_t *>(data);
              CGCHECK(cg_elements_read(cgnsFilePtr, base, zone, sect, idata, nullptr));
            }
            else {
              std::vector<cgsize_t> connect(element_nodes * num_to_get);
              CGCHECK(cg_elements_read(cgnsFilePtr, base, zone, sect, connect.data(), nullptr));
              if (field.get_type() == Ioss::Field::INT32) {
                auto * idata = reinterpret_cast<int *>(data);
                size_t i     = 0;
                for (auto node : connect) {
                  idata[i++] = node;
                }
              }
              else {
                auto * idata = reinterpret_cast<int64_t *>(data);
                size_t i     = 0;
                for (auto node : connect) {
                  idata[i++] = node;
                }
              }
            }
          }

          if (field.get_name() == "connectivity") {
            // Now need to map block-local node connectivity to global nodes...
            const auto &block_map = m_blockLocalNodeMap[zone];
            if (field.get_type() == Ioss::Field::INT32) {
              int *idata = static_cast<int *>(data);
              for (size_t i = 0; i < element_nodes * num_to_get; i++) {
                idata[i] = block_map[idata[i] - 1] + 1;
              }
            }
            else {
              int64_t *idata = static_cast<int64_t *>(data);
              for (size_t i = 0; i < element_nodes * num_to_get; i++) {
                idata[i] = block_map[idata[i] - 1] + 1;
              }
            }
          }
        }
        else if (field.get_name() == "ids") {
          // Map the local ids in this element block
          // (eb_offset+1...eb_offset+1+my_element_count) to global element ids.
          size_t eb_offset_plus_one = eb->get_offset() + 1;
          if (field.get_type() == Ioss::Field::INT64) {
            int64_t *idata = static_cast<int64_t *>(data);
            std::iota(idata, idata + my_element_count, eb_offset_plus_one);
          }
          else {
            assert(field.get_type() == Ioss::Field::INT32);
            int *idata = static_cast<int *>(data);
            std::iota(idata, idata + my_element_count, eb_offset_plus_one);
          }
        }
        else if (field.get_name() == "implicit_ids") {
          size_t eb_offset_plus_one = eb->get_offset() + 1;
          if (field.get_type() == Ioss::Field::INT64) {
            int64_t *idata = static_cast<int64_t *>(data);
            std::iota(idata, idata + my_element_count, eb_offset_plus_one);
          }
          else {
            assert(field.get_type() == Ioss::Field::INT32);
            int *idata = static_cast<int *>(data);
            std::iota(idata, idata + my_element_count, eb_offset_plus_one);
          }
        }
        else {
          num_to_get = Ioss::Utils::field_warning(eb, field, "input");
        }
      }
      else if (role == Ioss::Field::TRANSIENT) {
        // Locate the FlowSolution node corresponding to the correct state/step/time
        // TODO: do this at read_meta_data() and store...
        int step = get_region()->get_current_state();
        int solution_index =
            Utils::find_solution_index(cgnsFilePtr, base, zone, step, CG_CellCenter);

        double * rdata        = static_cast<double *>(data);
        cgsize_t range_min[1] = {1};
        cgsize_t range_max[1] = {my_element_count};

        auto var_type   = field.transformed_storage();
        int  comp_count = var_type->component_count();
        if (comp_count == 1) {
          CGCHECK(cg_field_read(cgnsFilePtr, base, zone, solution_index, field.get_name().c_str(),
                                CG_RealDouble, range_min, range_max, rdata));
        }
        else {
          std::vector<double> cgns_data(my_element_count);
          for (int i = 0; i < comp_count; i++) {
            char        field_suffix_separator = get_field_separator();
            std::string var_name =
                var_type->label_name(field.get_name(), i + 1, field_suffix_separator);

            CGCHECK(cg_field_read(cgnsFilePtr, base, zone, solution_index, var_name.c_str(),
                                  CG_RealDouble, range_min, range_max, cgns_data.data()));
            for (cgsize_t j = 0; j < my_element_count; j++) {
              rdata[comp_count * j + i] = cgns_data[j];
            }
          }
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(eb, field, "output");
      }
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    Ioss::Field::RoleType role = field.get_role();
    int                   base = sb->get_property("base").get_int();
    int                   zone = sb->get_property("zone").get_int();

    cgsize_t num_to_get = field.verify(data_size);

    cgsize_t rmin[3] = {0, 0, 0};
    cgsize_t rmax[3] = {0, 0, 0};

    bool cell_field = Utils::is_cell_field(field);
    if (cell_field) {
      assert(num_to_get == sb->get_property("cell_count").get_int());
      if (num_to_get > 0) {
        rmin[0] = sb->get_property("offset_i").get_int() + 1;
        rmin[1] = sb->get_property("offset_j").get_int() + 1;
        rmin[2] = sb->get_property("offset_k").get_int() + 1;

        rmax[0] = rmin[0] + sb->get_property("ni").get_int() - 1;
        rmax[1] = rmin[1] + sb->get_property("nj").get_int() - 1;
        rmax[2] = rmin[2] + sb->get_property("nk").get_int() - 1;
      }
    }
    else {
      // cell nodal field.
      assert(num_to_get == sb->get_property("node_count").get_int());
      if (num_to_get > 0) {
        rmin[0] = sb->get_property("offset_i").get_int() + 1;
        rmin[1] = sb->get_property("offset_j").get_int() + 1;
        rmin[2] = sb->get_property("offset_k").get_int() + 1;

        rmax[0] = rmin[0] + sb->get_property("ni").get_int();
        rmax[1] = rmin[1] + sb->get_property("nj").get_int();
        rmax[2] = rmin[2] + sb->get_property("nk").get_int();
      }
    }

    assert(num_to_get ==
           (rmax[0] - rmin[0] + 1) * (rmax[1] - rmin[1] + 1) * (rmax[2] - rmin[2] + 1));
    double *rdata = static_cast<double *>(data);

    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "mesh_model_coordinates_x") {
        CGCHECK(cg_coord_read(cgnsFilePtr, base, zone, "CoordinateX", CG_RealDouble, rmin, rmax,
                              rdata));
      }

      else if (field.get_name() == "mesh_model_coordinates_y") {
        CGCHECK(cg_coord_read(cgnsFilePtr, base, zone, "CoordinateY", CG_RealDouble, rmin, rmax,
                              rdata));
      }

      else if (field.get_name() == "mesh_model_coordinates_z") {
        CGCHECK(cg_coord_read(cgnsFilePtr, base, zone, "CoordinateZ", CG_RealDouble, rmin, rmax,
                              rdata));
      }

      else if (field.get_name() == "mesh_model_coordinates") {
        char basename[33];
        int  cell_dimension = 0;
        int  phys_dimension = 0;
        CGCHECK(cg_base_read(cgnsFilePtr, base, basename, &cell_dimension, &phys_dimension));

        // Data required by upper classes store x0, y0, z0, ... xn,
        // yn, zn. Data stored in cgns file is x0, ..., xn, y0,
        // ..., yn, z0, ..., zn so we have to allocate some scratch
        // memory to read in the data and then map into supplied
        // 'data'

        std::vector<double> coord(num_to_get);

        // ========================================================================
        // Repetitive code for each coordinate direction; use a lambda to consolidate...
        auto coord_lambda = [this, base, zone, &coord, rmin, rmax, phys_dimension, num_to_get,
                             &rdata](const char *ord_name, int ordinate) {
          CGCHECK(cg_coord_read(cgnsFilePtr, base, zone, ord_name, CG_RealDouble, rmin, rmax,
                                TOPTR(coord)));

          // Map to global coordinate position...
          for (cgsize_t i = 0; i < num_to_get; i++) {
            rdata[phys_dimension * i + ordinate] = coord[i];
          }
        };
        // End of lambda...
        // ========================================================================

        coord_lambda("CoordinateX", 0);

        if (phys_dimension >= 2) {
          coord_lambda("CoordinateY", 1);
        }

        if (phys_dimension == 3) {
          coord_lambda("CoordinateZ", 2);
        }
      }
      else if (field.get_name() == "cell_node_ids") {
        if (field.get_type() == Ioss::Field::INT64) {
          int64_t *idata = static_cast<int64_t *>(data);
          sb->get_cell_node_ids(idata, true);
        }
        else {
          assert(field.get_type() == Ioss::Field::INT32);
          int *idata = static_cast<int *>(data);
          sb->get_cell_node_ids(idata, true);
        }
      }
      else if (field.get_name() == "cell_ids") {
        if (field.get_type() == Ioss::Field::INT64) {
          int64_t *idata = static_cast<int64_t *>(data);
          sb->get_cell_ids(idata, true);
        }
        else {
          assert(field.get_type() == Ioss::Field::INT32);
          int *idata = static_cast<int *>(data);
          sb->get_cell_ids(idata, true);
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(sb, field, "input");
      }
    }
    else if (role == Ioss::Field::TRANSIENT) {
      auto var_type               = field.transformed_storage();
      int  comp_count             = var_type->component_count();
      char field_suffix_separator = get_field_separator();

      int sol_index = 0;
      int step      = get_region()->get_current_state();
      if (cell_field) {
        sol_index = Utils::find_solution_index(cgnsFilePtr, base, zone, step, CG_CellCenter);
      }
      else {
        sol_index = Utils::find_solution_index(cgnsFilePtr, base, zone, step, CG_Vertex);
      }

      if (comp_count == 1) {
        CGCHECK(cg_field_read(cgnsFilePtr, base, zone, sol_index, field.get_name().c_str(),
                              CG_RealDouble, rmin, rmax, rdata));
      }
      else {
        std::vector<double> cgns_data(num_to_get);
        for (int i = 0; i < comp_count; i++) {
          std::string var_name =
              var_type->label_name(field.get_name(), i + 1, field_suffix_separator);
          CGCHECK(cg_field_read(cgnsFilePtr, base, zone, sol_index, var_name.c_str(), CG_RealDouble,
                                rmin, rmax, cgns_data.data()));
          for (cgsize_t j = 0; j < num_to_get; j++) {
            rdata[comp_count * j + i] = cgns_data[j];
          }
        }
      }
    }
    else {
      num_to_get = Ioss::Utils::field_warning(sb, field, "input");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(ns, field, "input");
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(es, field, "input");
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(fs, field, "input");
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(es, field, "input");
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    int base = sb->get_property("base").get_int();
    int zone = sb->get_property("zone").get_int();
    int sect = sb->get_property("section").get_int();

    ssize_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {
      int64_t entity_count = sb->entity_count();
      if (num_to_get != entity_count) {
        std::ostringstream errmsg;
        errmsg << "ERROR: Partial field input not yet implemented for side blocks";
        IOSS_ERROR(errmsg);
      }
    }

    Ioss::Field::RoleType role = field.get_role();
    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "element_side_raw" || field.get_name() == "element_side") {

        // TODO(gdsjaar): ? Possibly rewrite using cgi_read_int_data so can skip reading element
        // connectivity
        int                   nodes_per_face = sb->topology()->number_nodes();
        std::vector<cgsize_t> elements(nodes_per_face * num_to_get); // Not needed, but can't skip

        // We get:
        // *  num_to_get parent elements,
        // *  num_to_get zeros (other parent element for face, but on boundary so 0)
        // *  num_to_get face_on_element
        // *  num_to_get zeros (face on other parent element)
        std::vector<cgsize_t> parent(4 * num_to_get);

        CGCHECK(cg_elements_read(cgnsFilePtr, base, zone, sect, TOPTR(elements), TOPTR(parent)));

        size_t offset = m_zoneOffset[zone - 1];
        if (field.get_type() == Ioss::Field::INT32) {
          int *  idata = reinterpret_cast<int *>(data);
          size_t j     = 0;
          for (ssize_t i = 0; i < num_to_get; i++) {
            idata[j++] = parent[num_to_get * 0 + i] + offset; // Element
            idata[j++] = parent[num_to_get * 2 + i];
            assert(parent[num_to_get * 1 + i] == 0);
            assert(parent[num_to_get * 3 + i] == 0);
          }
          // Adjust face numbers to IOSS convention instead of CGNS convention...
          Utils::map_cgns_face_to_ioss(sb->parent_element_topology(), num_to_get, idata);
        }
        else {
          int64_t *idata = reinterpret_cast<int64_t *>(data);
          size_t   j     = 0;
          for (ssize_t i = 0; i < num_to_get; i++) {
            idata[j++] = parent[num_to_get * 0 + i] + offset; // Element
            idata[j++] = parent[num_to_get * 2 + i];
            assert(parent[num_to_get * 1 + i] == 0);
            assert(parent[num_to_get * 3 + i] == 0);
          }
          // Adjust face numbers to IOSS convention instead of CGNS convention...
          Utils::map_cgns_face_to_ioss(sb->parent_element_topology(), num_to_get, idata);
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(sb, field, "input");
      }
    }
    else {
      num_to_get = Ioss::Utils::field_warning(sb, field, "input");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::get_field_internal(const Ioss::SideSet *fs, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(fs, field, "input");
  }
  int64_t DatabaseIO::get_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(cs, field, "input");
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::Region *region, const Ioss::Field &field,
                                         void * /*data*/, size_t /*data_size*/) const
  {
    return Ioss::Utils::field_warning(region, field, "output");
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::StructuredBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    Ioss::Field::RoleType role = field.get_role();
    int                   base = sb->get_property("base").get_int();
    int                   zone = sb->get_property("zone").get_int();

    cgsize_t num_to_get = field.verify(data_size);

    if (role == Ioss::Field::MESH) {
      bool cell_field = Utils::is_cell_field(field);

      if (cell_field) {
        assert(num_to_get == sb->get_property("cell_count").get_int());
      }

      double *rdata = static_cast<double *>(data);

      int crd_idx = 0;
      if (field.get_name() == "mesh_model_coordinates_x") {
        assert(!cell_field);
        CGCHECK(
            cg_coord_write(cgnsFilePtr, base, zone, CG_RealDouble, "CoordinateX", rdata, &crd_idx));
      }

      else if (field.get_name() == "mesh_model_coordinates_y") {
        assert(!cell_field);
        CGCHECK(
            cg_coord_write(cgnsFilePtr, base, zone, CG_RealDouble, "CoordinateY", rdata, &crd_idx));
      }

      else if (field.get_name() == "mesh_model_coordinates_z") {
        assert(!cell_field);
        CGCHECK(
            cg_coord_write(cgnsFilePtr, base, zone, CG_RealDouble, "CoordinateZ", rdata, &crd_idx));
      }

      else if (field.get_name() == "mesh_model_coordinates") {
        assert(!cell_field);
        int phys_dimension = get_region()->get_property("spatial_dimension").get_int();

        // Data required by upper classes store x0, y0, z0, ... xn,
        // yn, zn. Data stored in cgns file is x0, ..., xn, y0,
        // ..., yn, z0, ..., zn so we have to allocate some scratch
        // memory to read in the data and then map into supplied
        // 'data'
        std::vector<double> coord(num_to_get);

        // ========================================================================
        // Repetitive code for each coordinate direction; use a lambda to consolidate...
        auto coord_lambda = [this, &coord, num_to_get, phys_dimension, &rdata, base,
                             zone](const char *ord_name, int ordinate) {
          int crd_index = 0;

          // Map to global coordinate position...
          for (cgsize_t i = 0; i < num_to_get; i++) {
            coord[i] = rdata[phys_dimension * i + ordinate];
          }

          CGCHECK(cg_coord_write(cgnsFilePtr, base, zone, CG_RealDouble, ord_name, TOPTR(coord),
                                 &crd_index));
        };
        // End of lambda...
        // ========================================================================

        coord_lambda("CoordinateX", 0);

        if (phys_dimension >= 2) {
          coord_lambda("CoordinateY", 1);
        }
        if (phys_dimension == 3) {
          coord_lambda("CoordinateZ", 2);
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(sb, field, "output");
      }
    }
    else if (role == Ioss::Field::TRANSIENT) {
      double *rdata                  = static_cast<double *>(data);
      int     cgns_field             = 0;
      auto    var_type               = field.transformed_storage();
      int     comp_count             = var_type->component_count();
      char    field_suffix_separator = get_field_separator();
      if (comp_count == 1) {
        CGCHECK(cg_field_write(cgnsFilePtr, base, zone, m_currentCellCenterSolutionIndex,
                               CG_RealDouble, field.get_name().c_str(), rdata, &cgns_field));
        Utils::set_field_index(field, cgns_field, CG_CellCenter);
      }
      else {
        std::vector<double> cgns_data(num_to_get);
        for (int i = 0; i < comp_count; i++) {
          for (cgsize_t j = 0; j < num_to_get; j++) {
            cgns_data[j] = rdata[comp_count * j + i];
          }
          std::string var_name =
              var_type->label_name(field.get_name(), i + 1, field_suffix_separator);

          CGCHECK(cg_field_write(cgnsFilePtr, base, zone, m_currentCellCenterSolutionIndex,
                                 CG_RealDouble, var_name.c_str(), cgns_data.data(), &cgns_field));
          if (i == 0) {
            Utils::set_field_index(field, cgns_field, CG_CellCenter);
          }
        }
      }
    }
    else {
      num_to_get = Ioss::Utils::field_warning(sb, field, "output");
    }

    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::ElementBlock *eb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    size_t num_to_get = field.verify(data_size);
    if (num_to_get > 0) {

      Ioss::Field::RoleType role = field.get_role();

      if (role == Ioss::Field::MESH) {
        // Handle the MESH fields required for a CGNS file model.
        // (The 'genesis' portion)
        if (field.get_name() == "connectivity") {
          // This blocks zone has not been defined.
          // Get the "node block" for this element block...
          int element_nodes = eb->topology()->number_nodes();
          assert(field.raw_storage()->component_count() == element_nodes);

          Ioss::MapContainer nodes;
          nodes.reserve(element_nodes * num_to_get + 1);
          nodes.push_back(1); // Non-one-to-one map

          if (field.get_type() == Ioss::Field::INT32) {
            auto *idata = reinterpret_cast<int *>(data);
            for (size_t i = 0; i < element_nodes * num_to_get; i++) {
              nodes.push_back(idata[i]);
            }
          }
          else {
            auto *idata = reinterpret_cast<int64_t *>(data);
            for (size_t i = 0; i < element_nodes * num_to_get; i++) {
              nodes.push_back(idata[i]);
            }
          }
          Ioss::Utils::uniquify(nodes, true);
	  assert(nodes[0] == 1);

          // Now, we have the node count and cell count so we can create a zone...
          int      base    = 1;
          int      zone    = 0;
          cgsize_t size[3] = {0, 0, 0};
          size[1]          = eb->entity_count();
          size[0]          = nodes.size() - 1;

          CGCHECK(
              cg_zone_write(cgnsFilePtr, base, eb->name().c_str(), size, CG_Unstructured, &zone));
          eb->property_update("zone", zone);
          eb->property_update("id", zone);
          eb->property_update("guid", zone);
          eb->property_update("section", 1);
          eb->property_update("base", base);

          // Now we have a valid zone so can update some data structures...
          m_zoneOffset[zone]                = m_zoneOffset[zone - 1] + size[1];
          m_globalToBlockLocalNodeMap[zone] = new Ioss::Map("element", "unknown", myProcessor);
          m_globalToBlockLocalNodeMap[zone]->map().swap(nodes);
          m_globalToBlockLocalNodeMap[zone]->build_reverse_map();

          // Need to map global nodes to block-local node connectivity
          const auto &block_map = m_globalToBlockLocalNodeMap[zone];
          if (field.get_type() == Ioss::Field::INT32) {
            block_map->reverse_map_data((int *)data, field, num_to_get * element_nodes);
          }
          else {
            block_map->reverse_map_data((int64_t *)data, field, num_to_get * element_nodes);
          }

          if (eb->entity_count() > 0) {
            CG_ElementType_t type            = Utils::map_topology_to_cgns(eb->topology()->name());
            int              sect            = 0;
            int              field_byte_size = (field.get_type() == Ioss::Field::INT32) ? 32 : 64;
            if (field_byte_size == CG_SIZEOF_SIZE) {
              CGCHECK(cg_section_write(cgnsFilePtr, base, zone, "HexElements", type, 1, num_to_get,
                                       0, (cgsize_t *)data, &sect));
            }
            else {
              std::vector<cgsize_t> connect;
              connect.reserve(element_nodes * num_to_get);
              if (field.get_type() == Ioss::Field::INT32) {
                auto *idata = reinterpret_cast<int *>(data);
                for (size_t i = 0; i < element_nodes * num_to_get; i++) {
                  connect.push_back(idata[i]);
                }
              }
              else {
                auto *idata = reinterpret_cast<int64_t *>(data);
                for (size_t i = 0; i < element_nodes * num_to_get; i++) {
                  connect.push_back(idata[i]);
                }
              }
              CGCHECK(cg_section_write(cgnsFilePtr, base, zone, "HexElements", type, 1, num_to_get,
                                       0, connect.data(), &sect));
            }
            m_bcOffset[zone] += num_to_get;
            eb->property_update("section", sect);
          }
        }
        else {
          num_to_get = Ioss::Utils::field_warning(eb, field, "output");
        }
      }
      else if (role == Ioss::Field::TRANSIENT) {
        int     base                   = eb->get_property("base").get_int();
        int     zone                   = eb->get_property("zone").get_int();
        double *rdata                  = static_cast<double *>(data);
        int     cgns_field             = 0;
        auto    var_type               = field.transformed_storage();
        int     comp_count             = var_type->component_count();
        char    field_suffix_separator = get_field_separator();
        if (comp_count == 1) {
          CGCHECK(cg_field_write(cgnsFilePtr, base, zone, m_currentCellCenterSolutionIndex,
                                 CG_RealDouble, field.get_name().c_str(), rdata, &cgns_field));
          Utils::set_field_index(field, cgns_field, CG_CellCenter);
        }
        else {
          std::vector<double> cgns_data(num_to_get);
          for (int i = 0; i < comp_count; i++) {
            for (size_t j = 0; j < num_to_get; j++) {
              cgns_data[j] = rdata[comp_count * j + i];
            }
            std::string var_name =
                var_type->label_name(field.get_name(), i + 1, field_suffix_separator);

            CGCHECK(cg_field_write(cgnsFilePtr, base, zone, m_currentCellCenterSolutionIndex,
                                   CG_RealDouble, var_name.c_str(), cgns_data.data(), &cgns_field));
            if (i == 0) {
              Utils::set_field_index(field, cgns_field, CG_CellCenter);
            }
          }
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(eb, field, "output");
      }
    }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::FaceBlock *fb, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(fb, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeBlock *eb, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(eb, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::NodeBlock *nb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    // Instead of outputting a global nodeblock's worth of data,
    // the data is output a "zone" at a time.
    // The m_globalToBlockLocalNodeMap[zone] map is used (Ioss::Map pointer)
    // This map is built during the output of block connectivity,
    // so for cgns unstructured mesh, we need to output ElementBlock connectivity
    // prior to outputting nodal coordinates.
    for (const auto &z : m_globalToBlockLocalNodeMap) {
      if (z.second == nullptr) {
        std::ostringstream errmsg;
        errmsg << "ERROR: CGNS: The globalToBlockLocalNodeMap is not defined, so nodal fields "
                  "cannot be output.";
        IOSS_ERROR(errmsg);
      }
    }

    Ioss::Field::RoleType role       = field.get_role();
    int                   base       = 1;
    cgsize_t              num_to_get = field.verify(data_size);

    if (role == Ioss::Field::MESH) {
      if (field.get_name() == "mesh_model_coordinates" ||
          field.get_name() == "mesh_model_coordinates_x" ||
          field.get_name() == "mesh_model_coordinates_y" ||
          field.get_name() == "mesh_model_coordinates_z") {
        double *rdata = static_cast<double *>(data);

        if (field.get_name() == "mesh_model_coordinates") {
          int spatial_dim = nb->get_property("component_degree").get_int();
          for (const auto &block : m_globalToBlockLocalNodeMap) {
            auto zone = block.first;
            // NOTE: 'block_map' has one more entry than node_count.  First entry is for something
            // else.
            //       'block_map' is 1-based.
            const auto &        block_map = block.second;
            std::vector<double> x(block_map->map().size() - 1);
            std::vector<double> y(block_map->map().size() - 1);
            std::vector<double> z(block_map->map().size() - 1);

            for (size_t i = 0; i < block_map->map().size() - 1; i++) {
              auto global = block_map->map()[i + 1] - 1;
              x[i]        = rdata[global * spatial_dim + 0];
              if (spatial_dim > 1) {
                y[i] = rdata[global * spatial_dim + 1];
              }
              if (spatial_dim > 2) {
                z[i] = rdata[global * spatial_dim + 2];
              }
            }

            // Create the zone
            // Output this zones coordinates...
            int crd_idx = 0;
            CGCHECK(cg_coord_write(cgnsFilePtr, base, zone, CG_RealDouble, "CoordinateX", TOPTR(x),
                                   &crd_idx));

            if (spatial_dim > 1) {
              CGCHECK(cg_coord_write(cgnsFilePtr, base, zone, CG_RealDouble, "CoordinateY",
                                     TOPTR(y), &crd_idx));
            }

            if (spatial_dim > 2) {
              CGCHECK(cg_coord_write(cgnsFilePtr, base, zone, CG_RealDouble, "CoordinateZ",
                                     TOPTR(z), &crd_idx));
            }
          }
        }
        else {
          // Outputting only a single coordinate value...
          for (const auto &block : m_globalToBlockLocalNodeMap) {
            auto zone = block.first;
            // NOTE: 'block_map' has one more entry than node_count.  First entry is for something
            // else.
            //       'block_map' is 1-based.
            const auto &        block_map = block.second;
            std::vector<double> xyz(block_map->map().size() - 1);

            for (size_t i = 0; i < block_map->map().size() - 1; i++) {
              auto global = block_map->map()[i + 1] - 1;
              xyz[i]      = rdata[global];
            }

            std::string cgns_name = "Invalid";
            if (field.get_name() == "mesh_model_coordinates_x") {
              cgns_name = "CoordinateX";
            }
            else if (field.get_name() == "mesh_model_coordinates_y") {
              cgns_name = "CoordinateY";
            }
            else if (field.get_name() == "mesh_model_coordinates_z") {
              cgns_name = "CoordinateZ";
            }
            // Create the zone
            // Output this zones coordinates...
            int crd_idx = 0;
            CGCHECK(cg_coord_write(cgnsFilePtr, base, zone, CG_RealDouble, cgns_name.c_str(),
                                   TOPTR(xyz), &crd_idx));
          }
        }
      }
      else {
        num_to_get = Ioss::Utils::field_warning(nb, field, "output");
      }
    }
    else if (role == Ioss::Field::TRANSIENT) {
      double *rdata      = static_cast<double *>(data);
      int     cgns_field = 0;

      for (const auto &block : m_globalToBlockLocalNodeMap) {
        auto zone = block.first;
        // NOTE: 'block_map' has one more entry than node_count.
        // First entry is for something else.  'block_map' is
        // 1-based.
        const auto &        block_map = block.second;
        std::vector<double> blk_data(block_map->map().size() - 1);

        auto var_type   = field.transformed_storage();
        int  comp_count = var_type->component_count();

        if (comp_count == 1) {
          for (size_t j = 0; j < block_map->map().size() - 1; j++) {
            auto global = block_map->map()[j + 1] - 1;
            blk_data[j] = rdata[global];
          }
          CGCHECK(cg_field_write(cgnsFilePtr, base, zone, m_currentVertexSolutionIndex,
                                 CG_RealDouble, field.get_name().c_str(), blk_data.data(),
                                 &cgns_field));
          Utils::set_field_index(field, cgns_field, CG_Vertex);
        }
        else {
          char field_suffix_separator = get_field_separator();

          for (int i = 0; i < comp_count; i++) {
            for (size_t j = 0; j < block_map->map().size() - 1; j++) {
              auto global = block_map->map()[j + 1] - 1;
              blk_data[j] = rdata[comp_count * global + i];
            }
            std::string var_name =
                var_type->label_name(field.get_name(), i + 1, field_suffix_separator);
            CGCHECK(cg_field_write(cgnsFilePtr, base, zone, m_currentVertexSolutionIndex,
                                   CG_RealDouble, var_name.c_str(), blk_data.data(), &cgns_field));
            if (i == 0) {
              Utils::set_field_index(field, cgns_field, CG_Vertex);
            }
          }
        }
      }
    }
    else {
      num_to_get = Ioss::Utils::field_warning(nb, field, "output");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::NodeSet *ns, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(ns, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::EdgeSet *es, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(es, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::FaceSet *fs, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(fs, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::ElementSet *es, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(es, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::SideBlock *sb, const Ioss::Field &field,
                                         void *data, size_t data_size) const
  {
    const Ioss::EntityBlock *parent_block = sb->parent_block();
    if (parent_block == nullptr) {
      std::ostringstream errmsg;
      errmsg << "ERROR: CGNS: SideBlock " << sb->name()
             << " does not have a parent-block specified.  This is required for CGNS output.";
      IOSS_ERROR(errmsg);
    }

    int     base       = parent_block->get_property("base").get_int();
    int     zone       = parent_block->get_property("zone").get_int();
    ssize_t num_to_get = field.verify(data_size);

    Ioss::Field::RoleType role = field.get_role();

    if (role == Ioss::Field::MESH) {
      // Handle the MESH fields required for a CGNS file model.
      // (The 'genesis' portion)
      if (field.get_name() == "element_side") {
        // Get name from parent sideset...
        auto &name = sb->owner()->name();

        CG_ElementType_t type = Utils::map_topology_to_cgns(sb->topology()->name());
        int              sect = 0;

        cgsize_t cg_start = m_bcOffset[zone] + 1;
        cgsize_t cg_end   = m_bcOffset[zone] + num_to_get;
        m_bcOffset[zone] += num_to_get;

        // NOTE: Currently not writing the "ElementConnectivity" data for the
        //       boundary condition.  It isn't used in the read and don't have
        //       the data so would have to generate it.  This may cause problems
        //       with codes that use the downstream data if they base the BC off
        //       of the nodes instead of the element/side info.
        CGCHECK(cg_section_partial_write(cgnsFilePtr, base, zone, name.c_str(), type, cg_start,
                                         cg_end, 0, &sect));

        sb->property_update("section", sect);

        size_t                offset = m_zoneOffset[zone - 1];
        std::vector<cgsize_t> parent(4 * num_to_get);

        if (field.get_type() == Ioss::Field::INT32) {
          int *  idata = reinterpret_cast<int *>(data);
          size_t j     = 0;
          for (ssize_t i = 0; i < num_to_get; i++) {
            parent[num_to_get * 0 + i] = idata[j++] - offset; // Element
            parent[num_to_get * 2 + i] = idata[j++];
          }
          // Adjust face numbers to IOSS convention instead of CGNS convention...
          Utils::map_ioss_face_to_cgns(sb->parent_element_topology(), num_to_get, parent);
        }
        else {
          int64_t *idata = reinterpret_cast<int64_t *>(data);
          size_t   j     = 0;
          for (ssize_t i = 0; i < num_to_get; i++) {
            parent[num_to_get * 0 + i] = idata[j++] - offset; // Element
            parent[num_to_get * 2 + i] = idata[j++];
          }
          // Adjust face numbers to IOSS convention instead of CGNS convention...
          Utils::map_ioss_face_to_cgns(sb->parent_element_topology(), num_to_get, parent);
        }

        CGCHECK(cg_parent_data_write(cgnsFilePtr, base, zone, sect, TOPTR(parent)));
        return num_to_get;
      }

      num_to_get = Ioss::Utils::field_warning(sb, field, "output");
    }
    else {
      num_to_get = Ioss::Utils::field_warning(sb, field, "output");
    }
    return num_to_get;
  }

  int64_t DatabaseIO::put_field_internal(const Ioss::SideSet *ss, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(ss, field, "output");
  }
  int64_t DatabaseIO::put_field_internal(const Ioss::CommSet *cs, const Ioss::Field &field,
                                         void * /* data */, size_t /* data_size */) const
  {
    return Ioss::Utils::field_warning(cs, field, "output");
  }

  void DatabaseIO::write_results_meta_data() {}

  unsigned DatabaseIO::entity_field_support() const
  {
    return Ioss::NODEBLOCK | Ioss::ELEMENTBLOCK | Ioss::STRUCTUREDBLOCK | Ioss::NODESET |
           Ioss::SIDESET | Ioss::REGION;
  }

} // namespace Iocgns
