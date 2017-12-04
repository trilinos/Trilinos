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

#include <Ioss_Bar2.h>
#include <Ioss_Bar3.h>
#include <Ioss_Hex20.h>
#include <Ioss_Hex27.h>
#include <Ioss_Hex8.h>
#include <Ioss_Node.h>
#include <Ioss_Pyramid13.h>
#include <Ioss_Pyramid14.h>
#include <Ioss_Pyramid5.h>
#include <Ioss_Quad4.h>
#include <Ioss_Quad8.h>
#include <Ioss_Quad9.h>
#include <Ioss_StructuredBlock.h>
#include <Ioss_Tet10.h>
#include <Ioss_Tet4.h>
#include <Ioss_Tri3.h>
#include <Ioss_Tri4.h>
#include <Ioss_Tri6.h>
#include <Ioss_Unknown.h>
#include <Ioss_Wedge15.h>
#include <Ioss_Wedge18.h>
#include <Ioss_Wedge6.h>
#include <set>

#include <cgns/Iocgns_Utils.h>

#include <cgnsconfig.h>
#if CG_BUILD_PARALLEL
#include <pcgnslib.h>
#else
#include <cgnslib.h>
#endif

#define CGERR(funcall)                                                                             \
  if ((funcall) != CG_OK) {                                                                        \
    Iocgns::Utils::cgns_error(file_ptr, __FILE__, __func__, __LINE__, -1);                         \
  }

namespace {
  struct Range
  {
    Range(int a, int b) : m_beg(a < b ? a : b), m_end(a < b ? b : a), m_reversed(b < a) {}

    int  m_beg;
    int  m_end;
    bool m_reversed;
  };

  bool overlaps(const Range &a, const Range &b) { return a.m_beg <= b.m_end && b.m_beg <= a.m_end; }

  bool bc_overlaps(const Ioss::StructuredBlock *block, const Ioss::BoundaryCondition &bc)
  {
    int ordinal[3];
    ordinal[0] = block->get_property("ni").get_int();
    ordinal[1] = block->get_property("nj").get_int();
    ordinal[2] = block->get_property("nk").get_int();

    if (ordinal[0] == 0 && ordinal[1] == 0 && ordinal[2] == 0) {
      return false;
    }

    int offset[3];
    offset[0] = block->get_property("offset_i").get_int();
    offset[1] = block->get_property("offset_j").get_int();
    offset[2] = block->get_property("offset_k").get_int();

    // Note that block range is nodes and m_ordinal[] is cells, so need to add 1 to range.
    Range z_i(1 + offset[0], ordinal[0] + offset[0] + 1);
    Range z_j(1 + offset[1], ordinal[1] + offset[1] + 1);
    Range z_k(1 + offset[2], ordinal[2] + offset[2] + 1);

    Range gc_i(bc.m_rangeBeg[0], bc.m_rangeEnd[0]);
    Range gc_j(bc.m_rangeBeg[1], bc.m_rangeEnd[1]);
    Range gc_k(bc.m_rangeBeg[2], bc.m_rangeEnd[2]);

    return overlaps(z_i, gc_i) && overlaps(z_j, gc_j) && overlaps(z_k, gc_k);
  }

  Range subset_range(const Range &a, const Range &b)
  {
    Range ret(std::max(a.m_beg, b.m_beg), std::min(a.m_end, b.m_end));
    ret.m_reversed = a.m_reversed || b.m_reversed;
    return ret;
  }

  void bc_subset_range(const Ioss::StructuredBlock *block, Ioss::BoundaryCondition &bc)
  {
    if (bc_overlaps(block, bc)) {
      int ordinal[3];
      ordinal[0] = block->get_property("ni").get_int();
      ordinal[1] = block->get_property("nj").get_int();
      ordinal[2] = block->get_property("nk").get_int();

      int offset[3];
      offset[0] = block->get_property("offset_i").get_int();
      offset[1] = block->get_property("offset_j").get_int();
      offset[2] = block->get_property("offset_k").get_int();

      // NOTE: Updates the range in bc

      // Note that block range is nodes and m_ordinal[] is cells, so need to add 1 to range.
      Range z_i(1 + offset[0], ordinal[0] + offset[0] + 1);
      Range z_j(1 + offset[1], ordinal[1] + offset[1] + 1);
      Range z_k(1 + offset[2], ordinal[2] + offset[2] + 1);

      Range gc_i(bc.m_rangeBeg[0], bc.m_rangeEnd[0]);
      Range gc_j(bc.m_rangeBeg[1], bc.m_rangeEnd[1]);
      Range gc_k(bc.m_rangeBeg[2], bc.m_rangeEnd[2]);

      Range gc_ii = subset_range(z_i, gc_i);
      Range gc_jj = subset_range(z_j, gc_j);
      Range gc_kk = subset_range(z_k, gc_k);

      bc.m_rangeBeg[0] = gc_ii.m_reversed ? gc_ii.m_end : gc_ii.m_beg;
      bc.m_rangeEnd[0] = gc_ii.m_reversed ? gc_ii.m_beg : gc_ii.m_end;
      bc.m_rangeBeg[1] = gc_jj.m_reversed ? gc_jj.m_end : gc_jj.m_beg;
      bc.m_rangeEnd[1] = gc_jj.m_reversed ? gc_jj.m_beg : gc_jj.m_end;
      bc.m_rangeBeg[2] = gc_kk.m_reversed ? gc_kk.m_end : gc_kk.m_beg;
      bc.m_rangeEnd[2] = gc_kk.m_reversed ? gc_kk.m_beg : gc_kk.m_end;
    }
    else {
      bc.m_rangeBeg = {{0, 0, 0}};
      bc.m_rangeEnd = {{0, 0, 0}};
    }
  }

  int extract_trailing_int(char *name)
  {
    // 'name' consists of an arbitray number of characters followed by
    // zero or more digits.  Return the integer value of the contiguous
    // set of trailing digits.
    // Example: Name42 returns 42;  Name_52or_perhaps_3_43 returns 43.

    size_t len   = std::strlen(name);
    int    nstep = 0;
    int    mul   = 1;
    for (size_t d = len; d > 0; d--) {
      if (isdigit(name[d - 1])) {
        nstep += mul * (name[d - 1] - '0');
        mul *= 10;
      }
      else {
        break;
      }
    }
    return nstep;
  }
} // namespace

void Iocgns::Utils::cgns_error(int cgnsid, const char *file, const char *function, int lineno,
                               int processor)
{
  std::ostringstream errmsg;
  errmsg << "CGNS error '" << cg_get_error() << "' at line " << lineno << " in file '" << file
         << "' in function '" << function << "'";
  if (processor >= 0) {
    errmsg << " on processor " << processor;
  }
  errmsg << ". Please report to gdsjaar@sandia.gov if you need help.";
  if (cgnsid > 0) {
#if CG_BUILD_PARALLEL
    cgp_close(cgnsid);
#else
    cg_close(cgnsid);
#endif
  }
  IOSS_ERROR(errmsg);
}

CG_ZoneType_t Iocgns::Utils::check_zone_type(int cgnsFilePtr)
{
  // ========================================================================
  // Get the number of zones (element blocks) in the mesh...
  int base      = 1;
  int num_zones = 0;
  CGCHECKNP(cg_nzones(cgnsFilePtr, base, &num_zones));

  CG_ZoneType_t common_zone_type = CG_ZoneTypeNull;

  for (int zone = 1; zone <= num_zones; zone++) {
    CG_ZoneType_t zone_type;
    CGCHECKNP(cg_zone_type(cgnsFilePtr, base, zone, &zone_type));

    if (common_zone_type == CG_ZoneTypeNull) {
      common_zone_type = zone_type;
    }

    if (common_zone_type != zone_type) {
      std::ostringstream errmsg;
      errmsg << "ERROR: CGNS: Zone " << zone << " is not the same zone type as previous zones."
             << " This is currently not allowed or supported (hybrid mesh).";
      IOSS_ERROR(errmsg);
    }
  }
  return common_zone_type;
}

size_t Iocgns::Utils::index(const Ioss::Field &field) { return field.get_index() & 0xffffffff; }

void Iocgns::Utils::set_field_index(const Ioss::Field &field, size_t index,
                                    CG_GridLocation_t location)
{
  if (location == CG_CellCenter) {
    index |= CG_CELL_CENTER_FIELD_ID;
  }
  if (location == CG_Vertex) {
    index |= CG_VERTEX_FIELD_ID;
  }
  field.set_index(index);
}

bool Iocgns::Utils::is_cell_field(const Ioss::Field &field)
{
  size_t index = field.get_index();
  if (index & CG_VERTEX_FIELD_ID) {
    return false;
  }
  else if (index & CG_CELL_CENTER_FIELD_ID) {
    return true;
  }
  return !(field.get_name() == "mesh_model_coordinates" ||
           field.get_name() == "mesh_model_coordinates_x" ||
           field.get_name() == "mesh_model_coordinates_y" ||
           field.get_name() == "mesh_model_coordinates_z" ||
           field.get_name() == "cell_node_ids"); // Default to cell field...
}

void Iocgns::Utils::common_write_meta_data(int file_ptr, const Ioss::Region &region,
                                           std::vector<size_t> &zone_offset)
{
  // Make sure mesh is not hybrid...
  if (region.mesh_type() == Ioss::MeshType::HYBRID) {
    std::ostringstream errmsg;
    errmsg << "ERROR: CGNS: The mesh on region " << region.name() << " is of type 'hybrid'."
           << " This is currently not allowed or supported.";
    IOSS_ERROR(errmsg);
  }

  int base           = 0;
  int phys_dimension = region.get_property("spatial_dimension").get_int();
  CGERR(cg_base_write(file_ptr, "Base", phys_dimension, phys_dimension, &base));

  CGERR(cg_goto(file_ptr, base, "end"));
  CGERR(cg_descriptor_write("Information", "IOSS: CGNS Writer version -1"));

  // Output the sidesets as Family_t nodes
  const auto &sidesets = region.get_sidesets();
  for (const auto &ss : sidesets) {
    int fam = 0;
    CGERR(cg_family_write(file_ptr, base, ss->name().c_str(), &fam));

    int         bc_index = 0;
    CG_BCType_t bocotype = CG_BCTypeNull;
    if (ss->property_exists("bc_type")) {
      bocotype = (CG_BCType_t)ss->get_property("bc_type").get_int();
    }

    int64_t id = fam;
    if (ss->property_exists("id")) {
      id = ss->get_property("id").get_int();
    }

    CGERR(cg_fambc_write(file_ptr, base, fam, "FamBC", bocotype, &bc_index));
    CGERR(cg_goto(file_ptr, base, "Family_t", fam, nullptr));
    CGERR(cg_descriptor_write("FamBC_TypeId", std::to_string(bocotype).c_str()));
    CGERR(cg_descriptor_write("FamBC_TypeName", BCTypeName[bocotype]));
    CGERR(cg_descriptor_write("FamBC_UserId", std::to_string(id).c_str()));
    CGERR(cg_descriptor_write("FamBC_UserName", ss->name().c_str()));
  }

  // NOTE: Element Block zone write is deferred to put_field_internal so can
  // generate the node count based on connectivity traversal...

  const auto &structured_blocks = region.get_structured_blocks();
  for (const auto &sb : structured_blocks) {
    int      zone    = 0;
    cgsize_t size[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    size[3]          = sb->get_property("ni_global").get_int();
    size[4]          = sb->get_property("nj_global").get_int();
    size[5]          = sb->get_property("nk_global").get_int();

    size[0] = size[3] + 1;
    size[1] = size[4] + 1;
    size[2] = size[5] + 1;

    CGERR(cg_zone_write(file_ptr, base, sb->name().c_str(), size, CG_Structured, &zone));
    sb->property_update("zone", zone);
    sb->property_update("base", base);

    assert(zone > 0);
    zone_offset[zone] = zone_offset[zone - 1] + sb->get_property("cell_count").get_int();

    // Add GridCoordinates Node...
    int grid_idx = 0;
    CGERR(cg_grid_write(file_ptr, base, zone, "GridCoordinates", &grid_idx));

    // Transfer boundary condition nodes...

    // The bc.m_ownerRange argument needs to be the union of the size on all processors
    // Instead of requiring that of the caller, do the union in this routine.
    // TODO: Calculate it outside of the loop...
    std::vector<cgsize_t> bc_range(sb->m_boundaryConditions.size() * 6);
    size_t idx = 0;
    for (const auto &bc : sb->m_boundaryConditions) {
      for (size_t i=0; i < 3; i++) {
	bc_range[idx++] = -bc.m_rangeBeg[i];
      }
      for (size_t i=0; i < 3; i++) {
	bc_range[idx++] = bc.m_rangeEnd[i];
      }
    }

    region.get_database()->util().global_array_minmax(bc_range, Ioss::ParallelUtils::DO_MAX);

    idx = 0;
    for (idx = 0; idx < bc_range.size(); idx += 6) {
      bc_range[idx+0] = -bc_range[idx+0];
      bc_range[idx+1] = -bc_range[idx+1];
      bc_range[idx+2] = -bc_range[idx+2];
    }

    idx = 0;
    for (const auto &bc : sb->m_boundaryConditions) {
      int bc_idx = 0;
      CGERR(cg_boco_write(file_ptr, base, zone, bc.m_bcName.c_str(), CG_FamilySpecified,
                          CG_PointRange, 2, &bc_range[idx], &bc_idx));
      CGERR(cg_goto(file_ptr, base, sb->name().c_str(), 0, "ZoneBC_t", 1, bc.m_bcName.c_str(), 0,
                    "end"));
      CGERR(cg_famname_write(bc.m_bcName.c_str()));
      CGERR(cg_boco_gridlocation_write(file_ptr, base, zone, bc_idx, CG_Vertex));
      idx += 6;
    }

    // Transfer Zone Grid Connectivity...
    std::set<std::string> defined; // Will have entry if the zgc already output.
    for (const auto &zgc : sb->m_zoneConnectivity) {
      if (zgc.m_intraBlock) {
        continue;
      }
      // In parallel decomposition, can have multiple copies of a zgc; only output the first
      // instance
      if (defined.insert(zgc.m_connectionName).second) {
        int zgc_idx = 0;
        CGERR(cg_1to1_write(file_ptr, base, zone, zgc.m_connectionName.c_str(),
                            zgc.m_donorName.c_str(), &zgc.m_ownerRange[0], &zgc.m_donorRange[0],
                            &zgc.m_transform[0], &zgc_idx));
      }
    }
  }
}

std::string Iocgns::Utils::map_cgns_to_topology_type(CG_ElementType_t type)
{
  std::string topology = "unknown";
  switch (type) {
  case CG_NODE: topology = Ioss::Node::name; break;
  case CG_BAR_2: topology = Ioss::Bar2::name; break;
  case CG_BAR_3: topology = Ioss::Bar3::name; break;
  case CG_TRI_3: topology = Ioss::Tri3::name; break;
  case CG_TRI_6: topology = Ioss::Tri6::name; break;
  case CG_QUAD_4: topology = Ioss::Quad4::name; break;
  case CG_QUAD_8: topology = Ioss::Quad8::name; break;
  case CG_QUAD_9: topology = Ioss::Quad9::name; break;
  case CG_TETRA_4: topology = Ioss::Tet4::name; break;
  case CG_TETRA_10: topology = Ioss::Tet10::name; break;
  case CG_PYRA_5: topology = Ioss::Pyramid5::name; break;
  case CG_PYRA_13: topology = Ioss::Pyramid13::name; break;
  case CG_PYRA_14: topology = Ioss::Pyramid14::name; break;
  case CG_PENTA_6: topology = Ioss::Wedge6::name; break;
  case CG_PENTA_15: topology = Ioss::Wedge15::name; break;
  case CG_PENTA_18: topology = Ioss::Wedge18::name; break;
  case CG_HEXA_8: topology = Ioss::Hex8::name; break;
  case CG_HEXA_20: topology = Ioss::Hex20::name; break;
  case CG_HEXA_27: topology = Ioss::Hex27::name; break;
  default:
    std::cerr << "WARNING: Found topology of type " << cg_ElementTypeName(type)
              << " which is not currently supported.\n";
    topology = Ioss::Unknown::name;
  }
  return topology;
}

CG_ElementType_t Iocgns::Utils::map_topology_to_cgns(const std::string &name)
{
  CG_ElementType_t topo = CG_ElementTypeNull;
  if (name == Ioss::Node::name) {
    topo = CG_NODE;
  }
  else if (name == Ioss::Bar2::name) {
    topo = CG_BAR_2;
  }
  else if (name == Ioss::Bar3::name) {
    topo = CG_BAR_3;
  }
  else if (name == Ioss::Tri3::name) {
    topo = CG_TRI_3;
  }
  else if (name == Ioss::Tri6::name) {
    topo = CG_TRI_6;
  }
  else if (name == Ioss::Quad4::name) {
    topo = CG_QUAD_4;
  }
  else if (name == Ioss::Quad8::name) {
    topo = CG_QUAD_8;
  }
  else if (name == Ioss::Quad9::name) {
    topo = CG_QUAD_9;
  }
  else if (name == Ioss::Tet4::name) {
    topo = CG_TETRA_4;
  }
  else if (name == Ioss::Tet10::name) {
    topo = CG_TETRA_10;
  }
  else if (name == Ioss::Pyramid5::name) {
    topo = CG_PYRA_5;
  }
  else if (name == Ioss::Pyramid13::name) {
    topo = CG_PYRA_13;
  }
  else if (name == Ioss::Pyramid14::name) {
    topo = CG_PYRA_14;
  }
  else if (name == Ioss::Wedge6::name) {
    topo = CG_PENTA_6;
  }
  else if (name == Ioss::Wedge15::name) {
    topo = CG_PENTA_15;
  }
  else if (name == Ioss::Wedge18::name) {
    topo = CG_PENTA_18;
  }
  else if (name == Ioss::Hex8::name) {
    topo = CG_HEXA_8;
  }
  else if (name == Ioss::Hex20::name) {
    topo = CG_HEXA_20;
  }
  else if (name == Ioss::Hex27::name) {
    topo = CG_HEXA_27;
  }
  else {
    std::cerr << "WARNING: Found topology of type " << name
              << " which is not currently supported.\n";
  }
  return topo;
}

void Iocgns::Utils::write_flow_solution_metadata(int file_ptr, Ioss::Region *region, int state,
                                                 int *vertex_solution_index,
                                                 int *cell_center_solution_index)
{
  std::string c_name = "CellCenterSolutionAtStep";
  std::string v_name = "VertexSolutionAtStep";
  std::string step   = std::to_string(state);
  c_name += step;
  v_name += step;

  const auto &nblocks          = region->get_node_blocks();
  const auto &nblock           = nblocks[0];
  bool        has_nodal_fields = nblock->field_count(Ioss::Field::TRANSIENT) > 0;

  // Create a lambda to avoid code duplication for similar treatment
  // of structured blocks and element blocks.
  auto sol_lambda = [=](Ioss::GroupingEntity *block) {
    int base = block->get_property("base").get_int();
    int zone = block->get_property("zone").get_int();
    if (has_nodal_fields) {
      CGERR(cg_sol_write(file_ptr, base, zone, v_name.c_str(), CG_Vertex, vertex_solution_index));
      CGERR(
          cg_goto(file_ptr, base, "Zone_t", zone, "FlowSolution_t", *vertex_solution_index, "end"));
      CGERR(cg_gridlocation_write(CG_Vertex));
      CGERR(cg_descriptor_write("Step", step.c_str()));
    }
    if (block->field_count(Ioss::Field::TRANSIENT) > 0) {
      CGERR(cg_sol_write(file_ptr, base, zone, c_name.c_str(), CG_CellCenter,
                         cell_center_solution_index));
      CGERR(cg_goto(file_ptr, base, "Zone_t", zone, "FlowSolution_t", *cell_center_solution_index,
                    "end"));
      CGERR(cg_descriptor_write("Step", step.c_str()));
    }
  };

  // Use the lambda
  const auto &sblocks = region->get_structured_blocks();
  for (auto &block : sblocks) {
    sol_lambda(block);
  }
  // Use the lambda
  const auto &eblocks = region->get_element_blocks();
  for (auto &block : eblocks) {
    sol_lambda(block);
  }
}

int Iocgns::Utils::find_solution_index(int cgnsFilePtr, int base, int zone, int step,
                                       CG_GridLocation_t location)
{
  auto str_step = std::to_string(step);
  int  nsols    = 0;
  CGCHECKNP(cg_nsols(cgnsFilePtr, base, zone, &nsols));
  bool location_matches = false;
  for (int i = 0; i < nsols; i++) {
    CG_GridLocation_t db_location;
    char              db_name[33];
    CGCHECKNP(cg_sol_info(cgnsFilePtr, base, zone, i + 1, db_name, &db_location));
    if (location == db_location) {
      location_matches = true;
      // Check if steps match.
      // NOTE: Using non-standard "Descriptor_t" node in FlowSolution_t
      CGCHECKNP(cg_goto(cgnsFilePtr, base, "Zone_t", zone, "FlowSolution_t", i + 1, "end"));
      int descriptor_count = 0;
      CGCHECKNP(cg_ndescriptors(&descriptor_count));

      bool found_step_descriptor = false;
      for (int d = 0; d < descriptor_count; d++) {
        char *db_step = nullptr;
        char  name[33];
        CGCHECKNP(cg_descriptor_read(d + 1, name, &db_step));
        if (strcmp(name, "step") == 0) {
          found_step_descriptor = true;
          if (str_step == db_step) {
            cg_free(db_step);
            return i + 1;
          }

          cg_free(db_step);
          break; // Found "step" descriptor, but wasn't correct step...
        }
      }
      if (!found_step_descriptor) {
        // There was no Descriptor_t node with the name "step",
        // Try to decode the step from the FlowSolution_t name.
        int nstep = extract_trailing_int(db_name);
        if (nstep == step) {
          return i + 1;
        }
      }
    }
  }

  if (nsols == 1 && location_matches) {
    return 1;
  }
  std::cerr << "WARNING: CGNS: Could not find valid solution index for step " << step << ", zone "
            << zone << ", and location " << GridLocationName[location] << "\n";
  return 0;
}

void Iocgns::Utils::add_sidesets(int cgnsFilePtr, Ioss::DatabaseIO *db)
{
  int base         = 1;
  int num_families = 0;
  CGCHECKNP(cg_nfamilies(cgnsFilePtr, base, &num_families));

  for (int family = 1; family <= num_families; family++) {
    char        name[33];
    CG_BCType_t bocotype;
    int         num_bc  = 0;
    int         num_geo = 0;
    CGCHECKNP(cg_family_read(cgnsFilePtr, base, family, name, &num_bc, &num_geo));

#if IOSS_DEBUG_OUTPUT
    if (db->parallel_rank() == 0) {
      std::cerr << "Family " << family << " named " << name << " has " << num_bc << " BC, and "
                << num_geo << " geometry references\n";
    }
#endif
    if (num_bc > 0) {
      // Create a sideset...
      std::string ss_name(name); // Use name here before cg_fambc_read call overwrites it...

      CGCHECKNP(cg_fambc_read(cgnsFilePtr, base, family, 1, name, &bocotype));

      auto *ss = new Ioss::SideSet(db, ss_name);
      ss->property_add(Ioss::Property("id", family));
      ss->property_add(Ioss::Property("guid", db->util().generate_guid(family)));
      ss->property_add(Ioss::Property("bc_type", bocotype));
      db->get_region()->add(ss);
    }
  }
}

size_t Iocgns::Utils::resolve_nodes(Ioss::Region &region, int my_processor, bool is_parallel)
{
  // Each structured block has its own set of "cell_nodes"
  // At block boundaries, there are duplicate nodes which need to be resolved for the
  // unstructured mesh output.

  // We need to iterate all of the blocks and then each blocks zgc to determine
  // which nodes are shared between blocks. For all shared nodes, the node in the lowest
  // numbered zone is considered the "owner" and all other nodes are shared.

  // Create a vector of size which is the sum of the on-processor cell_nodes size for each block
  size_t num_total_cell_nodes = 0;
  auto & blocks               = region.get_structured_blocks();
  for (auto &block : blocks) {
    size_t node_count = block->get_property("node_count").get_int();
    num_total_cell_nodes += node_count;
  }

  ssize_t              ss_max = std::numeric_limits<ssize_t>::max();
  std::vector<ssize_t> cell_node_map(num_total_cell_nodes, ss_max);

  // Each cell_node location in the cell_node_map is currently initialized to ss_max.
  // Iterate each block and then each blocks intra-block (i.e., not
  // due to proc decomps) zgc instances and update cell_node_map
  // such that for each shared node, it points to the owner nodes
  // location.
  for (auto &block : blocks) {
    for (const auto &zgc : block->m_zoneConnectivity) {
      if (!zgc.m_intraBlock && zgc.m_isActive) { // Not due to processor decomposition.
        // NOTE: In parallel, the owner block should exist, but may not have
        // any cells on this processor.  We can access its global i,j,k, but
        // don't store or access any "bulk" data on it.
        auto owner_block = region.get_structured_block(zgc.m_donorName);
        assert(owner_block != nullptr);

        std::vector<int> i_range = zgc.get_range(1);
        std::vector<int> j_range = zgc.get_range(2);
        std::vector<int> k_range = zgc.get_range(3);
        for (auto &k : k_range) {
          for (auto &j : j_range) {
            for (auto &i : i_range) {
              Ioss::IJK_t index{{i, j, k}};
              Ioss::IJK_t owner = zgc.transform(index);

              // The nodes as 'index' and 'owner' are contiguous and
              // should refer to the same node. 'owner' should be
              // the owner (unless it is already owned by another
              // block)

              ssize_t global_offset = block->get_global_node_offset(index[0], index[1], index[2]);
              ssize_t owner_global_offset =
                  owner_block->get_global_node_offset(owner[0], owner[1], owner[2]);

              if (global_offset > owner_global_offset) {
                if (is_parallel && (zgc.m_donorProcessor != my_processor)) {
                  size_t block_local_offset =
                      block->get_block_local_node_offset(index[0], index[1], index[2]);
                  block->m_globalIdMap.emplace_back(block_local_offset, owner_global_offset + 1);
                }
                else if (!is_parallel || (zgc.m_ownerProcessor != my_processor)) {
                  size_t  local_offset = block->get_local_node_offset(index[0], index[1], index[2]);
                  ssize_t owner_local_offset =
                      owner_block->get_local_node_offset(owner[0], owner[1], owner[2]);

                  if (cell_node_map[local_offset] == ss_max) {
                    cell_node_map[local_offset] = owner_local_offset;
                  }
#if IOSS_DEBUG_OUTPUT
                  else {
                    if (cell_node_map[local_offset] != owner_local_offset) {
                      std::cerr << "DUPLICATE?: " << local_offset << " " << owner_local_offset
                                << " " << cell_node_map[local_offset] << " " << global_offset << " "
                                << owner_global_offset << "\n";
                    }
                  }
#endif
                }
              }
            }
          }
        }
      }
    }
  }

  // Now iterate cell_node_map.  If an entry == ss_max, then it is
  // an owned node and needs to have its index into the unstructed
  // mesh node map set; otherwise, the value points to the owner
  // node, so the index at this location should be set to the owner
  // nodes index.
  size_t index = 0;
  for (auto &node : cell_node_map) {
    if (node == ss_max || node < 0) {
      node = index++;
    }
    else {
      node = -node;
    }
  }

  for (auto &node : cell_node_map) {
    if (node < 0) {
      node = cell_node_map[-node];
    }
  }

  for (auto &block : blocks) {
    size_t node_count = block->get_property("node_count").get_int();
    block->m_blockLocalNodeIndex.resize(node_count);

    size_t beg = block->get_node_offset();
    size_t end = beg + node_count;
    for (size_t idx = beg, i = 0; idx < end; idx++) {
      block->m_blockLocalNodeIndex[i++] = cell_node_map[idx];
    }
  }
  return index;
}

void Iocgns::Utils::resolve_shared_nodes(Ioss::Region &region, int my_processor)
{
  // Determine which nodes are shared across processor boundaries.
  // Only need to check on block boundaries..

  // We need to iterate all of the blocks and then each blocks zgc to determine
  // which nodes are shared between processors. For all shared nodes, the node in the lowest
  // numbered zone is considered the "owner" and all other nodes are shared.

  // Create a vector of size which is the sum of the on-processor cell_nodes size for each block
  size_t num_total_cell_nodes = 0;
  auto & blocks               = region.get_structured_blocks();
  for (auto &block : blocks) {
    size_t node_count = block->get_property("node_count").get_int();
    num_total_cell_nodes += node_count;
  }

  ssize_t              ss_max = std::numeric_limits<ssize_t>::max();
  std::vector<ssize_t> cell_node_map(num_total_cell_nodes, ss_max);

  // Each cell_node location in the cell_node_map is currently
  // initialized to ss_max.  Iterate each block and then each blocks
  // intra-block (i.e., due to proc decomps) zgc instances and
  // update cell_node_map such that for each shared node, it points to
  // the owner nodes location.
  for (auto &owner_block : blocks) {
    auto owner_ids = owner_block->get_cell_node_ids(true);
    for (const auto &zgc : owner_block->m_zoneConnectivity) {
      if (zgc.m_isActive &&
          (zgc.m_donorProcessor != my_processor ||
           zgc.m_ownerProcessor != my_processor)) { // Due to processor decomposition.
        // NOTE: In parallel, the donor block should exist, but may not have
        // any cells on this processor.  We can access its global i,j,k, but
        // don't store or access any "bulk" data on it.
        auto donor_block = region.get_structured_block(zgc.m_donorName);
        assert(donor_block != nullptr);

        auto donor_ids = donor_block->get_cell_node_ids(true);

        std::vector<int> i_range = zgc.get_range(1);
        std::vector<int> j_range = zgc.get_range(2);
        std::vector<int> k_range = zgc.get_range(3);
        for (auto &k : k_range) {
          for (auto &j : j_range) {
            for (auto &i : i_range) {
              Ioss::IJK_t index{{i, j, k}};
              Ioss::IJK_t owner = zgc.transform(index);

              // The nodes as 'index' and 'owner' are contiguous and
              // should refer to the same node.

              if (my_processor == zgc.m_ownerProcessor) {
                ssize_t owner_offset =
                    owner_block->get_block_local_node_offset(index[0], index[1], index[2]);
                owner_block->m_sharedNode.emplace_back(owner_offset, zgc.m_donorProcessor);
              }
              else if (my_processor == zgc.m_donorProcessor) {
                ssize_t donor_offset =
                    donor_block->get_block_local_node_offset(owner[0], owner[1], owner[2]);
                donor_block->m_sharedNode.emplace_back(donor_offset, zgc.m_ownerProcessor);
              }
            }
          }
        }
      }
    }
#if 0 && IOSS_DEBUG_OUTPUT
    std::cerr << "P" << my_processor << ", Block " << owner_block->name()
              << " Shared Nodes: " << owner_block->m_sharedNode.size() << "\n";
#endif
  }
}

void Iocgns::Utils::add_structured_boundary_conditions(int                    cgnsFilePtr,
                                                       Ioss::StructuredBlock *block)
{
  int base = block->get_property("base").get_int();
  int zone = block->get_property("zone").get_int();
  int num_bcs;
  CGCHECKNP(cg_nbocos(cgnsFilePtr, base, zone, &num_bcs));

  cgsize_t range[6];
  for (int ibc = 0; ibc < num_bcs; ibc++) {
    char              boconame[33];
    CG_BCType_t       bocotype;
    CG_PointSetType_t ptset_type;
    cgsize_t          npnts;
    cgsize_t          NormalListSize;
    CG_DataType_t     NormalDataType;
    int               ndataset;

    // All we really want from this is 'boconame'
    CGCHECKNP(cg_boco_info(cgnsFilePtr, base, zone, ibc + 1, boconame, &bocotype, &ptset_type,
                           &npnts, nullptr, &NormalListSize, &NormalDataType, &ndataset));

    CGCHECKNP(cg_boco_read(cgnsFilePtr, base, zone, ibc + 1, range, nullptr));

    // There are some BC that are applied on an edge or a vertex;
    // Don't want those (yet?), so filter them out at this time...
    int same_count = (range[0] == range[3] ? 1 : 0) + (range[1] == range[4] ? 1 : 0) +
                     (range[2] == range[5] ? 1 : 0);
    if (same_count != 1) {
      std::cerr << "WARNING: CGNS: Skipping Boundary Condition '" << boconame << "' on block '"
                << block->name() << "'. It is applied to "
                << (same_count == 2 ? "an edge" : "a vertex")
                << ". This code only supports surfaces.\n";
      continue;
    }
    Ioss::SideSet *sset = block->get_database()->get_region()->get_sideset(boconame);
    if (sset == nullptr) {
      // Need to create a new sideset since didn't see this earlier.
      auto *db = block->get_database();
      sset     = new Ioss::SideSet(db, boconame);
      if (sset == nullptr) {
        std::ostringstream errmsg;
        errmsg << "ERROR: CGNS: Could not create sideset named '" << boconame << "' on block '"
               << block->name() << "'.\n";
        IOSS_ERROR(errmsg);
      }
      sset->property_add(Ioss::Property("id", ibc + 1)); // Not sure this is unique id...8
      sset->property_add(Ioss::Property("guid", db->util().generate_guid(ibc+1)));
      db->get_region()->add(sset);
    }

    if (sset != nullptr) {
      Ioss::IJK_t range_beg{{(int)std::min(range[0], range[3]), (int)std::min(range[1], range[4]),
                             (int)std::min(range[2], range[5])}};

      Ioss::IJK_t range_end{{(int)std::max(range[0], range[3]), (int)std::max(range[1], range[4]),
                             (int)std::max(range[2], range[5])}};

      // Determine overlap of surface with block (in parallel, a block may
      // be split among multiple processors and the block face this is applied
      // to may not exist on this decomposed block)
      auto        bc   = Ioss::BoundaryCondition(boconame, range_beg, range_end);
      std::string name = std::string(boconame) + "/" + block->name();

      bc_subset_range(block, bc);
      block->m_boundaryConditions.push_back(bc);
      auto sb = new Ioss::SideBlock(block->get_database(), name, Ioss::Quad4::name, Ioss::Hex8::name,
                                    block->m_boundaryConditions.back().get_face_count());
      sb->set_parent_block(block);
      sset->add(sb);
      sb->property_add(Ioss::Property("base", base));
      sb->property_add(Ioss::Property("zone", zone));
      sb->property_add(Ioss::Property("section", ibc + 1));
      sb->property_add(Ioss::Property("id", sset->get_property("id").get_int()));
      sb->property_add(Ioss::Property("guid", block->get_database()->util().generate_guid(sset->get_property("id").get_int())));

      // Set a property on the sideset specifying the boundary condition type (bocotype)
      // In CGNS, the bocotype is an enum; we store it as the integer value of the enum.
      if (sset->property_exists("bc_type")) {
        // Check that the 'bocotype' value matches the value of the property.
        auto old_bocotype = sset->get_property("bc_type").get_int();
        if (old_bocotype != bocotype && bocotype != CG_FamilySpecified) {
          IOSS_WARNING << "On sideset " << sset->name()
                       << ", the boundary condition type was previously set to " << old_bocotype
                       << " which does not match the current value of " << bocotype
                       << ". It will keep the old value.\n";
        }
      }
      else {
        sset->property_add(Ioss::Property("bc_type", (int)bocotype));
      }
    }
    else {
      std::ostringstream errmsg;
      errmsg << "ERROR: CGNS: StructuredBlock '" << block->name()
             << "' Did not find matching sideset with name '" << boconame << "'";
      IOSS_ERROR(errmsg);
    }
  }
}

void Iocgns::Utils::finalize_database(int cgnsFilePtr, const std::vector<double> &timesteps,
                                      Ioss::Region *region, int myProcessor)
{
  int base = 1;
  CGCHECK(cg_biter_write(cgnsFilePtr, base, "TimeIterValues", timesteps.size()));

  // Now write the timestep time values...
  CGCHECK(cg_goto(cgnsFilePtr, base, "BaseIterativeData_t", 1, "end"));
  cgsize_t dimtv[1] = {(cgsize_t)timesteps.size()};
  CGCHECK(cg_array_write("TimeValues", CG_RealDouble, 1, dimtv, timesteps.data()));

  // Output the ZoneIterativeData which maps a zones flow solutions to timesteps.
  // One per zone and the number of entries matches the number of timesteps...
  const auto &nblocks = region->get_node_blocks();
  auto &      nblock  = nblocks[0];

  bool has_nodal_fields = nblock->field_count(Ioss::Field::TRANSIENT) > 0;

  cgsize_t          dim[2] = {32, (cgsize_t)timesteps.size()};
  std::vector<char> names(32 * timesteps.size(), ' ');
  for (size_t state = 0; state < timesteps.size(); state++) {
    // This name is the "postfix" or common portion of all FlowSolution names...
    std::string name = "SolutionAtStep" + std::to_string(state + 1);
    std::strncpy(&names[state * 32], name.c_str(), 32);
    for (size_t i = name.size(); i < 32; i++) {
      names[state * 32 + i] = ' ';
    }
  }

  // Create a lambda to avoid code duplication for similar treatment
  // of structured blocks and element blocks.
  auto ziter = [=](Ioss::GroupingEntity *block) {
    int              zone = block->get_property("zone").get_int();
    std::vector<int> indices(timesteps.size());
    bool             has_cell_center_fields = block->field_count(Ioss::Field::TRANSIENT) > 0;
    if (has_cell_center_fields || has_nodal_fields) {
      CGCHECK(cg_ziter_write(cgnsFilePtr, base, zone, "ZoneIterativeData"));
      CGCHECK(cg_goto(cgnsFilePtr, base, "Zone_t", zone, "ZoneIterativeData_t", 1, "end"));
      CGCHECK(cg_array_write("FlowSolutionPointers", CG_Character, 2, dim, names.data()));

      if (has_nodal_fields) {
        int index     = 1;
        int increment = has_cell_center_fields ? 2 : 1;
        for (size_t state = 0; state < timesteps.size(); state++) {
          indices[state] = index;
          index += increment;
        }

        CGCHECK(cg_array_write("VertexSolutionIndices", CG_Integer, 1, &dim[1], indices.data()));
        CGCHECK(cg_descriptor_write("VertexPrefix", "Vertex"));
      }
      if (has_cell_center_fields) {
        int index     = has_nodal_fields ? 2 : 1;
        int increment = has_nodal_fields ? 2 : 1;
        for (size_t state = 0; state < timesteps.size(); state++) {
          indices[state] = index;
          index += increment;
        }

        CGCHECK(cg_array_write("CellCenterIndices", CG_Integer, 1, &dim[1], indices.data()));
        CGCHECK(cg_descriptor_write("CellCenterPrefix", "CellCenter"));
      }
    }
  };

  // Use the lambda...
  const auto &sblocks = region->get_structured_blocks();
  for (auto &block : sblocks) {
    ziter(block);
  }

  // Use the lambda...
  const auto &eblocks = region->get_element_blocks();
  for (auto &block : eblocks) {
    ziter(block);
  }
}

void Iocgns::Utils::add_transient_variables(int cgnsFilePtr, const std::vector<double> &timesteps,
                                            Ioss::Region *region, int myProcessor)
{
  // ==========================================
  // Add transient variables (if any) to all zones...
  // Create a lambda to avoid code duplication for similar treatment
  // of structured blocks and element blocks.

  // Assuming that the fields on all steps are the same, but can vary
  // from zone to zone.
  auto sol_iter = [=](Ioss::GroupingEntity *block) {
    int b = block->get_property("base").get_int();
    int z = block->get_property("zone").get_int();

    int sol_count = 0;
    CGCHECK(cg_nsols(cgnsFilePtr, b, z, &sol_count));
    int sol_per_step = sol_count / (int)timesteps.size();
    assert(sol_count % (int)timesteps.size() == 0);

    for (int sol = 1; sol <= sol_per_step; sol++) {
      char              solution_name[33];
      CG_GridLocation_t grid_loc;
      CGCHECK(cg_sol_info(cgnsFilePtr, b, z, sol, solution_name, &grid_loc));

      int field_count = 0;
      CGCHECK(cg_nfields(cgnsFilePtr, b, z, sol, &field_count));

      char **field_names = Ioss::Utils::get_name_array(field_count, 33);
      for (int field = 1; field <= field_count; field++) {
        CG_DataType_t data_type;
        char          field_name[33];
        CGCHECK(cg_field_info(cgnsFilePtr, b, z, sol, field, &data_type, field_name));
        std::strncpy(field_names[field - 1], field_name, 32);
      }

      // Convert raw field names into composite fields (a_x, a_y, a_z ==> 3D vector 'a')
      std::vector<Ioss::Field> fields;
      if (grid_loc == CG_CellCenter) {
        size_t entity_count = block->entity_count();
        Ioss::Utils::get_fields(entity_count, field_names, field_count, Ioss::Field::TRANSIENT, '_',
                                nullptr, fields);
        size_t index = 1;
        for (const auto &field : fields) {
          Utils::set_field_index(field, index, grid_loc);
          index += field.raw_storage()->component_count();
          block->field_add(field);
        }
      }
      else {
        assert(grid_loc == CG_Vertex);
        const Ioss::NodeBlock *cnb =
            (block->type() == Ioss::STRUCTUREDBLOCK)
                ? &(dynamic_cast<Ioss::StructuredBlock *>(block)->get_node_block())
                : region->get_node_blocks()[0];
        Ioss::NodeBlock *nb           = const_cast<Ioss::NodeBlock *>(cnb);
        size_t           entity_count = nb->entity_count();
        Ioss::Utils::get_fields(entity_count, field_names, field_count, Ioss::Field::TRANSIENT, '_',
                                nullptr, fields);
        size_t index = 1;
        for (const auto &field : fields) {
          Utils::set_field_index(field, index, grid_loc);
          index += field.raw_storage()->component_count();
          nb->field_add(field);
        }
      }

      Ioss::Utils::delete_name_array(field_names, field_count);
    }
  };
  // ==========================================

  if (!timesteps.empty()) {
    const auto &sblocks = region->get_structured_blocks();
    for (auto &block : sblocks) {
      sol_iter(block);
    }
    const auto &eblocks = region->get_element_blocks();
    for (auto &block : eblocks) {
      sol_iter(block);
    }
  }
}

int Iocgns::Utils::get_step_times(int cgnsFilePtr, std::vector<double> &timesteps,
                                  Ioss::Region *region, double timeScaleFactor, int myProcessor)
{
  int  base          = 1;
  int  num_timesteps = 0;
  char bitername[33];
  int  ierr = cg_biter_read(cgnsFilePtr, base, bitername, &num_timesteps);
  if (ierr == CG_NODE_NOT_FOUND) {
    return num_timesteps;
  }
  if (ierr == CG_ERROR) {
    Utils::cgns_error(cgnsFilePtr, __FILE__, __func__, __LINE__, myProcessor);
  }

  if (num_timesteps <= 0) {
    return num_timesteps;
  }

  // Read the timestep time values.
  CGCHECK(cg_goto(cgnsFilePtr, base, "BaseIterativeData_t", 1, "end"));
  std::vector<double> times(num_timesteps);
  CGCHECK(cg_array_read_as(1, CG_RealDouble, times.data()));

  timesteps.reserve(num_timesteps);
  for (int i = 0; i < num_timesteps; i++) {
    region->add_state(times[i] * timeScaleFactor);
    timesteps.push_back(times[i]);
  }
  return num_timesteps;
}
