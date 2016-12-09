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

#include <cgns/Iocgns_Utils.h>

#include <cgnsconfig.h>
#if CG_BUILD_PARALLEL
#include <pcgnslib.h>
#else
#include <cgnslib.h>
#endif

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
}

void Iocgns::Utils::cgns_error(int cgnsid, const char *file, const char *function, int lineno,
                               int processor)
{
  std::ostringstream errmsg;
  errmsg << "CGNS error '" << cg_get_error() << "' at line " << lineno << " in file '" << file
         << "' in function '" << function << "' on processor " << processor
         << ". Please report to gdsjaar@sandia.gov if you need help.";
  if (cgnsid > 0) {
#if CG_BUILD_PARALLEL
    cgp_close(cgnsid);
#else
    cg_close(cgnsid);
#endif
  }
  IOSS_ERROR(errmsg);
}

std::string Iocgns::Utils::map_cgns_to_topology_type(CG_ElementType_t type)
{
  std::string topology = "unknown";
  switch (type) {
  case CG_NODE: topology     = Ioss::Node::name; break;
  case CG_BAR_2: topology    = Ioss::Bar2::name; break;
  case CG_BAR_3: topology    = Ioss::Bar3::name; break;
  case CG_TRI_3: topology    = Ioss::Tri3::name; break;
  case CG_TRI_6: topology    = Ioss::Tri6::name; break;
  case CG_QUAD_4: topology   = Ioss::Quad4::name; break;
  case CG_QUAD_8: topology   = Ioss::Quad8::name; break;
  case CG_QUAD_9: topology   = Ioss::Quad9::name; break;
  case CG_TETRA_4: topology  = Ioss::Tet4::name; break;
  case CG_TETRA_10: topology = Ioss::Tet10::name; break;
  case CG_PYRA_5: topology   = Ioss::Pyramid5::name; break;
  case CG_PYRA_13: topology  = Ioss::Pyramid13::name; break;
  case CG_PYRA_14: topology  = Ioss::Pyramid14::name; break;
  case CG_PENTA_6: topology  = Ioss::Wedge6::name; break;
  case CG_PENTA_15: topology = Ioss::Wedge15::name; break;
  case CG_PENTA_18: topology = Ioss::Wedge18::name; break;
  case CG_HEXA_8: topology   = Ioss::Hex8::name; break;
  case CG_HEXA_20: topology  = Ioss::Hex20::name; break;
  case CG_HEXA_27: topology  = Ioss::Hex27::name; break;
  default:
    std::cerr << "WARNING: Found topology of type " << cg_ElementTypeName(type)
              << " which is not currently supported.\n";
    topology = Ioss::Unknown::name;
  }
  return topology;
}

void Iocgns::Utils::add_sidesets(int cgnsFilePtr, Ioss::DatabaseIO *db)
{
  cgsize_t base         = 1;
  cgsize_t num_families = 0;
  cg_nfamilies(cgnsFilePtr, base, &num_families);
  for (cgsize_t family = 1; family <= num_families; family++) {
    char        name[33];
    CG_BCType_t bocotype;
    cgsize_t    num_bc  = 0;
    cgsize_t    num_geo = 0;
    cg_family_read(cgnsFilePtr, base, family, name, &num_bc, &num_geo);
#if defined(IOSS_DEBUG_OUTPUT)
    std::cout << "Family " << family << " named " << name << " has " << num_bc << " BC, and "
              << num_geo << " geometry references\n";
#endif
    if (num_bc > 0) {
      // Create a sideset...
      std::string ss_name(name); // Use name here before cg_fambc_read call overwrites it...

      cg_fambc_read(cgnsFilePtr, base, family, 1, name, &bocotype);
      auto *ss = new Ioss::SideSet(db, ss_name);
      ss->property_add(Ioss::Property("id", family));
      ss->property_add(Ioss::Property("bc_type", bocotype));
      db->get_region()->add(ss);
    }
  }
}

size_t Iocgns::Utils::resolve_nodes(Ioss::Region &region, int my_processor)
{
  // Each structured block has its own set of "cell_nodes"
  // At block boundaries, there are duplicate nodes which need to be resolved for the
  // unstructured mesh output.

  // We need to iterate all of the blocks and then each blocks zgc to determine
  // which nodes are shared between blocks. For all shared nodes, the node in the lowest
  // numbered zone is considered the "owner" and all other nodes are shared.

  // Create a vector of size which is the sum of the on-processor cell_nodes size for each block
  auto &blocks = region.get_structured_blocks();

  ssize_t ss_max               = std::numeric_limits<ssize_t>::max();
  size_t  num_total_cell_nodes = 0;
  for (auto &block : blocks) {
    size_t node_count = block->get_property("node_count").get_int();
    num_total_cell_nodes += node_count;
  }
  std::vector<ssize_t> cell_node_map(num_total_cell_nodes, ss_max);

  // Each cell_node location in the cell_node_map is currently initialized to ss_max.
  // Iterate each block and then each blocks intra-block (i.e., not
  // due to proc decomps) zgc instances and update cell_node_map
  // such that for each shared node, it points to the owner nodes
  // location.
  for (auto &block : blocks) {
    for (const auto &zgc : block->m_zoneConnectivity) {
      if (!zgc.m_intraBlock) { // Not due to processor decomposition.
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
                assert(zgc.m_donorProcessor != -1);
                if (zgc.m_donorProcessor != my_processor) {
                  size_t block_local_offset =
                      block->get_block_local_node_offset(index[0], index[1], index[2]);
                  block->m_globalIdMap.emplace_back(block_local_offset, owner_global_offset + 1);
                }
                else {
                  size_t  local_offset = block->get_local_node_offset(index[0], index[1], index[2]);
                  ssize_t owner_local_offset =
                      owner_block->get_local_node_offset(owner[0], owner[1], owner[2]);

                  if (cell_node_map[local_offset] == ss_max) {
                    cell_node_map[local_offset] = owner_local_offset;
                  }
#if defined(IOSS_DEBUG_OUTPUT)
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

void Iocgns::Utils::add_structured_boundary_conditions(int                    cgnsFilePtr,
                                                       Ioss::StructuredBlock *block)
{
  cgsize_t base = block->get_property("base").get_int();
  cgsize_t zone = block->get_property("zone").get_int();
  int      num_bcs;
  cg_nbocos(cgnsFilePtr, base, zone, &num_bcs);

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
    cg_boco_info(cgnsFilePtr, base, zone, ibc + 1, boconame, &bocotype, &ptset_type, &npnts,
                 nullptr, &NormalListSize, &NormalDataType, &ndataset);

    cg_boco_read(cgnsFilePtr, base, zone, ibc + 1, range, nullptr);

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
      sset->property_add(Ioss::Property("id", ibc + 1)); // Not sure this is unique id...
      db->get_region()->add(sset);
    }

    assert(sset != nullptr);
    Ioss::IJK_t range_beg{{std::min(range[0], range[3]), std::min(range[1], range[4]),
	  std::min(range[2], range[5])}};

    Ioss::IJK_t range_end{{std::max(range[0], range[3]), std::max(range[1], range[4]),
	  std::max(range[2], range[5])}};

    // Determine overlap of surface with block (in parallel, a block may
    // be split among multiple processors and the block face this is applied
    // to may not exist on this decomposed block)
    auto        bc   = Ioss::BoundaryCondition(boconame, range_beg, range_end);
    std::string name = std::string(boconame) + "/" + block->name();

    if (bc_overlaps(block, bc)) {
      bc_subset_range(block, bc);
      block->m_boundaryConditions.push_back(bc);
      auto sb = new Ioss::SideBlock(block->get_database(), name, "Quad4", "Hex8",
				    block->m_boundaryConditions.back().get_face_count());
      sb->set_parent_block(block);
      sset->add(sb);
      sb->property_add(Ioss::Property("base", base));
      sb->property_add(Ioss::Property("zone", zone));
      sb->property_add(Ioss::Property("section", ibc + 1));
    }
    else {
      Ioss::IJK_t zeros{{0, 0, 0}};
      auto        zero_bc = Ioss::BoundaryCondition(boconame, zeros, zeros);
      block->m_boundaryConditions.push_back(zero_bc);
      auto sb = new Ioss::SideBlock(block->get_database(), name, "Quad4", "Hex8", 0);
      sb->set_parent_block(block);
      sset->add(sb);
      sb->property_add(Ioss::Property("base", base));
      sb->property_add(Ioss::Property("zone", zone));
      sb->property_add(Ioss::Property("section", ibc + 1));
    }

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
}
