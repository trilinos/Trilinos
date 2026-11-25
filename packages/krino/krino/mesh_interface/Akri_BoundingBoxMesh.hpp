// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_BoundingBoxMesh_h
#define Akri_BoundingBoxMesh_h

#include <Akri_MeshInterface.hpp>
#include <Akri_BoundingBox.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <array>
#include <memory>

namespace krino {

class CartesianCoordinateMapping
{
public:
  CartesianCoordinateMapping(const size_t nx, const size_t ny, const size_t nz,
      const BoundingBox_T<double,3> & bbox)
  : m_N{nx, ny, nz},
    m_bbox(bbox)
     {}
  
  std::array<double, 3> get_node_coordinates(const int spatial_dim, const std::array<size_t, 3> & ijk) const
  {
    std::array<double, 3> node_coords{};
    const stk::math::Vector3d & min = m_bbox.get_min();
    const stk::math::Vector3d & max = m_bbox.get_max();
    for (int i=0; i < spatial_dim; ++i)
      node_coords[i] = min[i] + (max[i]-min[i])*ijk[i]/m_N[i];
    return node_coords;
  }

  std::array<double, 3> get_triangle_lattice_node_coordinates(const bool flattenBoundaries, const std::array<size_t, 2> & ij) const
  {
    std::array<double, 3> node_coords{};
    const stk::math::Vector3d & min = m_bbox.get_min();
    const stk::math::Vector3d & max = m_bbox.get_max();
    const double offset = ij[1]%2==0 ? 0.0 : -0.5;
    node_coords[0] = min[0] + (max[0]-min[0])/m_N[0] * (ij[0]+offset);
    if (flattenBoundaries)
        node_coords[0] = std::min(max[0], std::max(min[0], node_coords[0]));
    node_coords[1] = min[1] + (max[1]-min[1])/m_N[1] * ij[1];
    return node_coords;
  }
  
  std::array<double, 3> get_BCC_node_coordinates(const int spatial_dim, const bool flattenBoundaries, const std::array<size_t, 3> & ijk, const std::array<int, 3> & dijk) const
  {
    std::array<double, 3> node_coords{};
    const stk::math::Vector3d & min = m_bbox.get_min();
    const stk::math::Vector3d & max = m_bbox.get_max();
    for (int i=0; i < spatial_dim; ++i)
    {
      node_coords[i] = min[i] + (max[i]-min[i])*(0.5+ijk[i]+dijk[i])/m_N[i];
      if (flattenBoundaries)
        node_coords[i] = std::min(max[i], std::max(min[i], node_coords[i]));
    }
    return node_coords;
  }
private:
  const std::array<size_t,3> m_N;
  const BoundingBox_T<double,3> m_bbox;
};

enum BoundingBoxMeshStructureType
{
  CUBIC_BOUNDING_BOX_MESH = 0,
  BCC_BOUNDING_BOX_MESH = 1,
  FLAT_WALLED_BCC_BOUNDING_BOX_MESH = 2,
  TRIANGULAR_LATTICE_BOUNDING_BOX_MESH = 3,
  FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH = 4
};

class BoundingBoxMesh : public MeshInterface {
public:
  typedef BoundingBox_T<double,3> BoundingBoxType;
  static std::array<size_t,3> get_node_x_y_z( stk::mesh::EntityId entity_id, const std::array<size_t,3> & N );

public:
  BoundingBoxMesh(stk::topology element_topology, stk::ParallelMachine comm);
  virtual void populate_mesh(const stk::mesh::BulkData::AutomaticAuraOption auto_aura_option = stk::mesh::BulkData::AUTO_AURA) override;
  virtual stk::mesh::MetaData & meta_data() override { STK_ThrowAssert( nullptr != m_meta.get() ) ; return *m_meta; }
  virtual const stk::mesh::MetaData & meta_data() const override { STK_ThrowAssert( nullptr != m_meta.get() ) ; return *m_meta; }
  virtual stk::mesh::BulkData & bulk_data() override { STK_ThrowAssert( nullptr != m_mesh.get() ) ; return *m_mesh; }
  virtual const stk::mesh::BulkData & bulk_data() const override { STK_ThrowAssert( nullptr != m_mesh.get() ) ; return *m_mesh; }

  void set_domain(const BoundingBoxType & mesh_bbox, const double mesh_size, const int pad_cells = 0);
  void create_domain_sides();
  const CartesianCoordinateMapping & get_coord_mapping() const { return *my_coord_mapping; }
  void get_node_x_y_z( stk::mesh::EntityId entity_id, size_t &ix , size_t &iy , size_t &iz ) const;
  void set_mesh_structure_type(BoundingBoxMeshStructureType type) { myMeshStructureType = type; }
  bool has_flat_boundaries() const { return CUBIC_BOUNDING_BOX_MESH == myMeshStructureType || FLAT_WALLED_BCC_BOUNDING_BOX_MESH == myMeshStructureType || FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH == myMeshStructureType; }
  // Number the nodes and elements in reverse (id -> maximum id + 1 - id), so that every ordering by
  // id is reversed while the mesh, and its decomposition across processors, stay the same.  For
  // testing that results do not depend on the ids.  Set before populate_mesh().  The static
  // get_node_x_y_z takes the unreversed ids.
  void set_reverse_entity_ids(const bool reverse) { myReverseEntityIds = reverse; }
  // Decompose the mesh across processors differently: in slabs along x (rows in y for a triangular
  // lattice) instead of by coordinate bisection, given to the processors in reverse order, so that
  // the processors' parts, and which processor owns each shared node, differ.  For testing that
  // results do not depend on the decomposition.  Set before populate_mesh().
  void set_alternate_processor_decomposition(const bool alternate) { myAlternateProcessorDecomposition = alternate; }
  // Whether populate_mesh puts the sides of the flat walls in the sidesets Xminus, Xplus, ... (by default
  // it does).  Without, the walls' sides are only in the exposed-boundary part, so that the walls are in
  // one surface; the sidesets are still declared, empty.  Set before populate_mesh().
  void set_create_domain_sidesets(const bool createDomainSidesets) { myCreateDomainSidesets = createDomainSidesets; }
private:
  // The id in the mesh of a node or element with the given (unreversed) id, and vice versa.
  stk::mesh::EntityId mesh_node_id(const stk::mesh::EntityId nodeId) const { return myReverseEntityIds ? myMaxNodeId + 1 - nodeId : nodeId; }
  stk::mesh::EntityId mesh_element_id(const stk::mesh::EntityId elemId) const { return myReverseEntityIds ? myMaxElementId + 1 - elemId : elemId; }
  static BoundingBoxMeshStructureType default_structure_type_for_topology(const stk::topology elementTopology);
  void declare_domain_side_parts(const stk::mesh::Part & blockPart);
  void require_has_flat_boundaries() const;
  void populate_2D_triangular_lattice_based_mesh();
  void populate_BCC_mesh();
  void populate_cell_based_mesh();
  void setup_cell_node( size_t ix , size_t iy , size_t iz );
  void setup_BCC_node( size_t ix , size_t iy , size_t iz, int dx, int dy, int dz );
  void build_face_tets( size_t cell_id, size_t ix , size_t iy , size_t iz, int iface, const std::vector<stk::mesh::EntityId> & cell_node_ids,
      const std::vector<std::array<size_t,3>> & cell_node_lattice_indices );
  bool is_quadratic() const { return m_element_topology != m_element_topology.base(); }
  // The mid-edge nodes of quadratic elements: each vertex has an integer index on a doubled lattice,
  // from 0 to maxVertexLatticeIndex in each direction, and the mid-edge node of an edge gets the id
  // of the sum of its vertices' indices (its midpoint on a lattice 4 times finer), after the
  // vertex ids.  No two edges of these meshes have the same midpoint, so the ids are unique, and
  // they are the same however the mesh is decomposed.
  void set_mid_edge_node_lattice(const std::array<size_t,3> & maxVertexLatticeIndex, const stk::mesh::EntityId maxVertexId);
  stk::mesh::EntityId get_mid_edge_node_id(const std::array<size_t,3> & vertex0LatticeIndex, const std::array<size_t,3> & vertex1LatticeIndex) const;
  // The element, with its mid-edge nodes if it is quadratic.
  void declare_element_with_mid_edge_nodes(const stk::mesh::EntityId elemId,
      const std::vector<stk::mesh::EntityId> & vertexIds,
      const std::vector<std::array<size_t,3>> & vertexLatticeIndices);
  void set_mid_edge_node_coordinates();
  stk::mesh::EntityId get_node_id(const size_t ix, const size_t iy, const size_t iz ) const
  {
    const size_t node_id_start = 1;
    return node_id_start + ix + ( m_nx + 1 ) * ( iy + ( m_ny + 1 ) * iz );
  }
  stk::mesh::EntityId get_BCC_node_id(const size_t ix, const size_t iy, const size_t iz, const int dx, const int dy, const int dz ) const
  {
    const size_t node_id_start = 1;
    const size_t num_nodes = ( m_nx + 1 ) * ( m_ny + 1 ) * ( m_nz + 1 );
    return node_id_start + num_nodes + (ix+1+dx) + ( m_nx + 2 ) * ( (iy+1+dy) + ( m_ny + 2 ) * (iz+1+dz) );
  }
  // Half-open rectangular range of structured cell indices [lo, hi) owned by a processor.
  struct CartesianCellRange
  {
    std::array<size_t,3> lo;
    std::array<size_t,3> hi;
    size_t size() const { return (hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]); }
  };
  // Recursive coordinate bisection of the structured cell grid. Deterministic and
  // communication-free: every processor computes the same decomposition, so any processor
  // can query the cell box owned by any other rank (needed to build node sharing).
  static CartesianCellRange recursive_coordinate_bisection(std::array<size_t,3> lo, std::array<size_t,3> hi,
      const int rank_begin, const int group_size, const int target_rank);
  CartesianCellRange determine_processor_cell_range(const int p_size, const int p_rank) const;
  void generate_node_to_processor_map(const int p_size,
      const int p_rank,
      const std::vector<std::array<int,3>> & cell_node_locations,
      std::unordered_map<stk::mesh::EntityId,
      std::vector<int>> & nodes_to_procs) const;
  void set_is_cell_edge_function_for_BCC_mesh() const;
  void set_is_cell_edge_function_for_cell_based_mesh() const;
private:
  stk::ParallelMachine myComm;
  std::shared_ptr<stk::mesh::MetaData> m_meta;
  std::unique_ptr<stk::mesh::BulkData> m_mesh;
  std::unique_ptr<CartesianCoordinateMapping> my_coord_mapping;
  stk::mesh::PartVector m_elem_parts;
  stk::mesh::PartVector m_node_parts;
  const stk::topology m_element_topology;
  BoundingBoxMeshStructureType myMeshStructureType{CUBIC_BOUNDING_BOX_MESH};
  size_t m_nx, m_ny, m_nz;
  std::array<size_t,3> myMidEdgeNodeLatticeSize{{0,0,0}};
  stk::mesh::EntityId myFirstMidEdgeNodeId{0};
  bool myReverseEntityIds{false};
  bool myAlternateProcessorDecomposition{false};
  stk::mesh::EntityId myMaxNodeId{0};
  stk::mesh::EntityId myMaxElementId{0};
  BoundingBoxType m_mesh_bbox;
  std::vector<stk::mesh::Part*> mySideParts;
  bool myCreateDomainSidesets{true};
};

} // namespace krino

#endif // Akri_BoundingBoxMesh_h
