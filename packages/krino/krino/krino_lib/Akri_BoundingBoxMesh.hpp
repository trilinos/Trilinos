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
  void get_node_coordinates(double * node_coords, const int spatial_dim, const std::array<size_t, 3> & ijk) const
  {
    const stk::math::Vector3d & min = m_bbox.get_min();
    const stk::math::Vector3d & max = m_bbox.get_max();
    for (int i=0; i < spatial_dim; ++i)
      node_coords[i] = min[i] + (max[i]-min[i])*ijk[i]/m_N[i];
  }
  void get_triangle_lattice_node_coordinates(double * node_coords, const bool flattenBoundaries, const std::array<size_t, 2> & ij) const
  {
    const stk::math::Vector3d & min = m_bbox.get_min();
    const stk::math::Vector3d & max = m_bbox.get_max();
    const double offset = ij[1]%2==0 ? 0.0 : -0.5;
    node_coords[0] = min[0] + (max[0]-min[0])/m_N[0] * (ij[0]+offset);
    if (flattenBoundaries)
        node_coords[0] = std::min(max[0], std::max(min[0], node_coords[0]));
    node_coords[1] = min[1] + (max[1]-min[1])/m_N[1] * ij[1];
  }
  void get_BCC_node_coordinates(double * node_coords, const int spatial_dim, const bool flattenBoundaries, const std::array<size_t, 3> & ijk, const std::array<int, 3> & dijk) const
  {
    const stk::math::Vector3d & min = m_bbox.get_min();
    const stk::math::Vector3d & max = m_bbox.get_max();
    for (int i=0; i < spatial_dim; ++i)
    {
      node_coords[i] = min[i] + (max[i]-min[i])*(0.5+ijk[i]+dijk[i])/m_N[i];
      if (flattenBoundaries)
        node_coords[i] = std::min(max[i], std::max(min[i], node_coords[i]));
    }
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
private:
  static BoundingBoxMeshStructureType default_structure_type_for_topology(const stk::topology elementTopology);
  void declare_domain_side_parts(const stk::mesh::Part & blockPart);
  void require_has_flat_boundaries() const;
  void populate_2D_triangular_lattice_based_mesh();
  void populate_BCC_mesh();
  void populate_cell_based_mesh();
  void setup_cell_node( size_t ix , size_t iy , size_t iz );
  void setup_BCC_node( size_t ix , size_t iy , size_t iz, int dx, int dy, int dz );
  void build_face_tets( size_t cell_id, size_t ix , size_t iy , size_t iz, int iface, const std::vector<stk::mesh::EntityId> & cell_node_ids );
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
  void get_cell_x_y_z( stk::mesh::EntityId cell_id, size_t &ix , size_t &iy , size_t &iz ) const;
  std::pair<size_t, size_t> determine_processor_cells(const int p_size, const int p_rank) const;
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
  BoundingBoxType m_mesh_bbox;
  std::vector<stk::mesh::Part*> mySideParts;
};

} // namespace krino

#endif // Akri_BoundingBoxMesh_h
