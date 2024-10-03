// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_AuxMetaData.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_DiagWriter.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_tools/mesh_tools/FixNodeSharingViaSearch.hpp>

namespace krino{

BoundingBoxMesh::BoundingBoxMesh(stk::topology element_topology, stk::ParallelMachine comm)
: myComm(comm),
  m_element_topology(element_topology), m_nx(0), m_ny(0), m_nz(0)
{
  STK_ThrowRequire(element_topology == stk::topology::TRIANGLE_3_2D ||
      element_topology == stk::topology::QUADRILATERAL_4_2D ||
      element_topology == stk::topology::TETRAHEDRON_4 ||
      element_topology == stk::topology::HEXAHEDRON_8);

  myMeshStructureType = default_structure_type_for_topology(element_topology);

  m_meta = stk::mesh::MeshBuilder().set_spatial_dimension(element_topology.dimension())
                                   .create_meta_data();

  AuxMetaData & aux_meta = AuxMetaData::create(*m_meta);
  stk::mesh::Part & block_part = m_meta->declare_part_with_topology( "block_1", element_topology );
  stk::io::put_io_part_attribute(block_part);
  m_elem_parts.push_back(&block_part);
  m_elem_parts.push_back(&aux_meta.active_part());
  m_node_parts.push_back(&aux_meta.active_part());

  declare_domain_side_parts(block_part);

  stk::mesh::Field<double> & coordsField = m_meta->declare_field<double>(stk::topology::NODE_RANK, "coordinates", 1);
  stk::mesh::put_field_on_mesh(coordsField, m_meta->universal_part(), m_meta->spatial_dimension(), nullptr);
}

BoundingBoxMeshStructureType BoundingBoxMesh::default_structure_type_for_topology(const stk::topology elementTopology)
{
  if (elementTopology == stk::topology::QUADRILATERAL_4_2D || elementTopology == stk::topology::HEXAHEDRON_8)
    return CUBIC_BOUNDING_BOX_MESH;
  else if (elementTopology == stk::topology::TETRAHEDRON_4)
    return FLAT_WALLED_BCC_BOUNDING_BOX_MESH;
  else if (elementTopology == stk::topology::TRIANGLE_3_2D)
    return FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH;
  STK_ThrowRequireMsg(false, "Unsupport topology " << elementTopology.name());
  return CUBIC_BOUNDING_BOX_MESH;
}

void
BoundingBoxMesh::set_domain(const BoundingBoxType & mesh_bbox, const double mesh_size, const int pad_size)
{
  m_mesh_bbox = mesh_bbox;
  const stk::math::Vector3d padding(0.5*mesh_size*pad_size, 0.5*mesh_size*pad_size, 0.5*mesh_size*pad_size);
  m_mesh_bbox = BoundingBoxType(mesh_bbox.get_min() - padding, mesh_bbox.get_max() + padding);

  const typename BoundingBoxType::VecType min = m_mesh_bbox.get_min();
  const typename BoundingBoxType::VecType max = m_mesh_bbox.get_max();
  const typename BoundingBoxType::VecType span = max-min;
  m_nx = (size_t) (0.5 + span[0] / mesh_size);
  const bool isTriangularLattice = (TRIANGULAR_LATTICE_BOUNDING_BOX_MESH == myMeshStructureType || FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH == myMeshStructureType);
  const double dy = isTriangularLattice ? (0.5*mesh_size*std::sqrt(3.0)) : mesh_size;
  m_ny = (size_t) (0.5 + span[1] / dy);
  m_nz = (m_element_topology.dimension() == 2) ? 1 : ((size_t) (0.5 + span[2] / mesh_size));

  if (m_element_topology.dimension() == 2)
  {
    krinolog << "Generated mesh with for domain "
      << "(min,max) = ((" << min[0] << "," << min[1] << "),("<< max[0] << "," << max[1] << ")"
      << ", (nx,ny) = (" << m_nx << "," << m_ny << ")" << stk::diag::dendl;
  }
  else
  {
    krinolog << "Generated mesh with for domain "
      << "(min,max) = ((" << min[0] << "," << min[1] << "," << min[2] << "),("<< max[0] << "," << max[1] << "," << max[2] << ")"
      << ", (nx,ny,nz) = (" << m_nx << "," << m_ny << "," << m_nz << ")" << stk::diag::dendl;
  }

}

void BoundingBoxMesh::set_is_cell_edge_function_for_cell_based_mesh() const
{
  const std::array<size_t,3> N = {m_nx, m_ny, m_nz};
  const int dim = m_meta->spatial_dimension();
  STK_ThrowRequireMsg(m_mesh, "Cannot call set_is_cell_edge_function() before BulkData is created.");
  const stk::mesh::BulkData & mesh = *m_mesh;

  auto is_cell_edge = [N,dim,&mesh](stk::mesh::Entity node0, stk::mesh::Entity node1)
  {
    const auto node0Indices = get_node_x_y_z(mesh.identifier(node0), N);
    const auto node1Indices = get_node_x_y_z(mesh.identifier(node1), N);
    if (dim == 2)
    {
      return node0Indices[0] == node1Indices[0] || node0Indices[1] == node1Indices[1];
    }
    const bool matchX = node0Indices[0] == node1Indices[0];
    const bool matchY = node0Indices[1] == node1Indices[1];
    const bool matchZ = node0Indices[2] == node1Indices[2];
    return (matchX && matchY) || (matchX && matchZ) || (matchY && matchZ);
  };

  AuxMetaData & aux_meta = AuxMetaData::get(*m_meta);
  aux_meta.set_is_cell_edge_function(is_cell_edge);
}

void BoundingBoxMesh::set_is_cell_edge_function_for_BCC_mesh() const
{
  auto is_cell_edge = [](stk::mesh::Entity node0, stk::mesh::Entity node1)
  {
    return true;
  };

  AuxMetaData & aux_meta = AuxMetaData::get(*m_meta);
  aux_meta.set_is_cell_edge_function(is_cell_edge);
}

void
BoundingBoxMesh::populate_mesh(const stk::mesh::BulkData::AutomaticAuraOption auto_aura_option)
{
  STK_ThrowRequireMsg(m_mesh_bbox.valid(), "Must call set_domain() before populate_mesh()");
  m_mesh = stk::mesh::MeshBuilder(myComm).set_aura_option(auto_aura_option).create(m_meta);
  if (CUBIC_BOUNDING_BOX_MESH == myMeshStructureType)
    populate_cell_based_mesh();
  else if (TRIANGULAR_LATTICE_BOUNDING_BOX_MESH == myMeshStructureType || FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH == myMeshStructureType)
    populate_2D_triangular_lattice_based_mesh();
  else if (BCC_BOUNDING_BOX_MESH == myMeshStructureType || FLAT_WALLED_BCC_BOUNDING_BOX_MESH == myMeshStructureType)
    populate_BCC_mesh();
  else
    STK_ThrowRequireMsg(false, "Unsupported or unrecognized mesh structure type " << myMeshStructureType);


  stk::mesh::create_exposed_block_boundary_sides(*m_mesh, m_meta->universal_part(), {&AuxMetaData::get(*m_meta).exposed_boundary_part()});

  if (has_flat_boundaries())
    create_domain_sides();
}

enum BCCNode { BCC_NODE=8, BCC_NODE_XMINUS=9, BCC_NODE_XPLUS=10, BCC_NODE_YMINUS=11, BCC_NODE_YPLUS=12, BCC_NODE_ZMINUS=13, BCC_NODE_ZPLUS=14 };

void
BoundingBoxMesh::build_face_tets( size_t cell_id, size_t ix , size_t iy , size_t iz, int iface, const std::vector<stk::mesh::EntityId> & cell_node_ids )
{
  std::vector<std::vector<std::vector<int>>> faceTets = {
      {{BCC_NODE, BCC_NODE_XMINUS, 0, 4},
       {BCC_NODE, BCC_NODE_XMINUS, 4, 7},
       {BCC_NODE, BCC_NODE_XMINUS, 7, 3},
       {BCC_NODE, BCC_NODE_XMINUS, 3, 0}},
      {{BCC_NODE, BCC_NODE_XPLUS, 1, 2},
       {BCC_NODE, BCC_NODE_XPLUS, 2, 6},
       {BCC_NODE, BCC_NODE_XPLUS, 6, 5},
       {BCC_NODE, BCC_NODE_XPLUS, 5, 1}},
      {{BCC_NODE, BCC_NODE_YMINUS, 0, 1},
       {BCC_NODE, BCC_NODE_YMINUS, 1, 5},
       {BCC_NODE, BCC_NODE_YMINUS, 5, 4},
       {BCC_NODE, BCC_NODE_YMINUS, 4, 0}},
      {{BCC_NODE, BCC_NODE_YPLUS, 2, 3},
       {BCC_NODE, BCC_NODE_YPLUS, 3, 7},
       {BCC_NODE, BCC_NODE_YPLUS, 7, 6},
       {BCC_NODE, BCC_NODE_YPLUS, 6, 2}},
      {{BCC_NODE, BCC_NODE_ZMINUS, 0, 3},
       {BCC_NODE, BCC_NODE_ZMINUS, 3, 2},
       {BCC_NODE, BCC_NODE_ZMINUS, 2, 1},
       {BCC_NODE, BCC_NODE_ZMINUS, 1, 0}},
      {{BCC_NODE, BCC_NODE_ZPLUS, 4, 5},
       {BCC_NODE, BCC_NODE_ZPLUS, 5, 6},
       {BCC_NODE, BCC_NODE_ZPLUS, 6, 7},
       {BCC_NODE, BCC_NODE_ZPLUS, 7, 4}}
  };

  const stk::mesh::EntityId elem_id_start = 1;
  const int num_elem_per_face = 4;
  const int num_faces_per_cell = 6;
  const int num_nodes_per_elem = 4;
  std::vector<stk::mesh::EntityId> elemNodes(num_nodes_per_elem);

  const auto & faceTet = faceTets[iface];
  for (int ielem = 0; ielem < num_elem_per_face; ++ielem)
  {
    stk::mesh::EntityId elem_id = num_elem_per_face*num_faces_per_cell*cell_id + num_elem_per_face*iface + ielem + elem_id_start; // Uses twice as many IDs as necessary
    for (int elem_node_index = 0; elem_node_index<num_nodes_per_elem; ++elem_node_index)
    {
      elemNodes[elem_node_index] = cell_node_ids[faceTet[ielem][elem_node_index]];
    }

    stk::mesh::declare_element( *m_mesh, m_elem_parts, elem_id, elemNodes );
  }
}

void
BoundingBoxMesh::setup_cell_node( size_t ix , size_t iy , size_t iz )
{
  stk::mesh::FieldBase const* coord_field = m_meta->coordinate_field();
  stk::mesh::EntityId nodeId = get_node_id(ix, iy, iz);
  stk::mesh::Entity const node = m_mesh->get_entity( stk::topology::NODE_RANK, nodeId );
  if (m_mesh->is_valid(node))
  {
    m_mesh->change_entity_parts(node, m_node_parts);

    double * coord_data = static_cast<double*>(stk::mesh::field_data( *coord_field, node ));
    my_coord_mapping->get_node_coordinates(coord_data, m_meta->spatial_dimension(), {ix, iy, iz});
  }
}

void
BoundingBoxMesh::setup_BCC_node( size_t ix , size_t iy , size_t iz, int dx, int dy, int dz )
{
  stk::mesh::FieldBase const* coord_field = m_meta->coordinate_field();
  stk::mesh::EntityId nodeId = get_BCC_node_id(ix, iy, iz, dx, dy, dz);
  stk::mesh::Entity const node = m_mesh->get_entity( stk::topology::NODE_RANK, nodeId );
  if (m_mesh->is_valid(node))
  {
    m_mesh->change_entity_parts(node, m_node_parts);

    const bool flattenBoundaries = FLAT_WALLED_BCC_BOUNDING_BOX_MESH == myMeshStructureType;
    double * coord_data = static_cast<double*>(stk::mesh::field_data( *coord_field, node ));
    my_coord_mapping->get_BCC_node_coordinates(coord_data, m_meta->spatial_dimension(), flattenBoundaries, {ix, iy, iz}, {dx, dy, dz});
  }
}

size_t get_triangle_lattice_node_id(size_t ix, size_t iy, size_t nx, size_t ny)
{
  return 1 + iy*(nx+1) + ix + ((iy%2 == 0) ? (iy/2) : ((iy-1)/2));
}

void fill_triangle_lattice_node_indices(const size_t ix, const size_t iy, std::array<std::array<size_t,2>,3> & triNodeIndices)
{
  if (iy%2 == 0)
  {
    if (ix%2 == 0)
      triNodeIndices = {{ {{ix/2,iy}}, {{ix/2+1,iy+1}}, {{ix/2,iy+1}} }};
    else
      triNodeIndices = {{ {{(ix-1)/2,iy}}, {{(ix-1)/2+1,iy}}, {{(ix-1)/2+1,iy+1}} }};
  }
  else
  {
    if (ix%2 == 0)
      triNodeIndices = {{ {{ix/2,iy}}, {{ix/2+1,iy}}, {{ix/2,iy+1}} }};
    else
      triNodeIndices = {{ {{(ix-1)/2+1,iy}}, {{(ix-1)/2+1,iy+1}}, {{(ix-1)/2,iy+1}} }};
  }
}

void
BoundingBoxMesh::populate_2D_triangular_lattice_based_mesh()
{ /* %TRACE[ON]% */ Trace trace__("krino::BoundingBoxMesh::populate_2D_triangular_lattice_based_mesh()"); /* %TRACE% */
  STK_ThrowRequire(m_mesh);

  const int p_size = m_mesh->parallel_size();
  const int p_rank = m_mesh->parallel_rank();
  stk::mesh::FieldBase const* coord_field = m_meta->coordinate_field();
  my_coord_mapping = std::make_unique<CartesianCoordinateMapping>(m_nx, m_ny, m_nz, m_mesh_bbox);

  STK_ThrowRequire(m_element_topology == stk::topology::TRIANGLE_3_2D);
  const bool flattenBoundaries = FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH == myMeshStructureType;

  const size_t NelemPerSlabInY = (2*m_nx+1);
  const size_t Nelem = m_ny*NelemPerSlabInY;
  const size_t beg_elem = ( Nelem * p_rank ) / p_size ;
  const size_t end_elem = ( Nelem * ( p_rank + 1 ) ) / p_size ;

  const size_t num_local_cells = end_elem-beg_elem;
  krinolog << "BoundingBoxMesh creating " << num_local_cells << " local cells." << stk::diag::dendl;

  std::vector<stk::mesh::EntityId> triNodes(3);
  std::array<std::array<size_t,2>,3> triNodeIndices;

  m_mesh->modification_begin();

  {
    for (size_t elemIndex=beg_elem; elemIndex!=end_elem; ++elemIndex)
    {
      const size_t iy = elemIndex / NelemPerSlabInY;
      const size_t ix = elemIndex - iy*NelemPerSlabInY;

      fill_triangle_lattice_node_indices(ix, iy, triNodeIndices);

      triNodes = {get_triangle_lattice_node_id(triNodeIndices[0][0], triNodeIndices[0][1], m_nx, m_ny),
          get_triangle_lattice_node_id(triNodeIndices[1][0], triNodeIndices[1][1], m_nx, m_ny),
          get_triangle_lattice_node_id(triNodeIndices[2][0], triNodeIndices[2][1], m_nx, m_ny)};

      stk::mesh::declare_element( *m_mesh, m_elem_parts, elemIndex+1, triNodes );

      for (int n=0; n<3; ++n)
      {
        stk::mesh::Entity const node = m_mesh->get_entity( stk::topology::NODE_RANK , triNodes[n] );
        m_mesh->change_entity_parts(node, m_node_parts);

        double * coord_data = static_cast<double*>(stk::mesh::field_data( *coord_field, node ));

        my_coord_mapping->get_triangle_lattice_node_coordinates(coord_data, flattenBoundaries, {{triNodeIndices[n][0], triNodeIndices[n][1]}});
      }
    }
    stk::tools::fix_node_sharing_via_search(*m_mesh);
  }
  m_mesh->modification_end();

  std::vector<size_t> counts;
  stk::mesh::count_entities( m_meta->locally_owned_part(), *m_mesh, counts );
  krinolog << "Generated mesh with " << counts[stk::topology::ELEM_RANK] << " local elements and " << counts[stk::topology::NODE_RANK] << " local nodes." << stk::diag::dendl;
}

void
BoundingBoxMesh::populate_BCC_mesh()
{ /* %TRACE[ON]% */ Trace trace__("krino::BoundingBoxMesh::populate_BCC_mesh()"); /* %TRACE% */
  STK_ThrowRequire(m_mesh);
  set_is_cell_edge_function_for_BCC_mesh();

  const int p_size = m_mesh->parallel_size();
  const int p_rank = m_mesh->parallel_rank();

  std::vector<std::array<int,3>> hex_cell_node_locations = { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1} };
  std::vector<std::pair<BCCNode,std::array<int,3>>> BCC_cell_neighbors = { {BCC_NODE,{0,0,0}},
      {BCC_NODE_XMINUS,{-1,0,0}}, {BCC_NODE_XPLUS,{+1,0,0}},
      {BCC_NODE_YMINUS,{0,-1,0}}, {BCC_NODE_YPLUS,{0,+1,0}},
      {BCC_NODE_ZMINUS,{0,0,-1}}, {BCC_NODE_ZPLUS,{0,0,+1}}};

  my_coord_mapping = std::make_unique<CartesianCoordinateMapping>(m_nx, m_ny, m_nz, m_mesh_bbox);

  std::pair<size_t,size_t> proc_cell_range = determine_processor_cells(p_size, p_rank);
  const size_t num_local_cells = proc_cell_range.second-proc_cell_range.first;
  krinolog << "BoundingBoxMesh creating " << num_local_cells << " local cells." << stk::diag::dendl;

  const int num_nodes_per_cell = 8;
  const int num_nodes_per_cell_plus_BCC_nodes = 15;
  std::vector<stk::mesh::EntityId> cell_node_ids(num_nodes_per_cell_plus_BCC_nodes);

  m_mesh->modification_begin();

  {
    size_t count = 0;
    for (size_t cell_id=proc_cell_range.first; cell_id!=proc_cell_range.second; ++cell_id)
    {
      size_t ix = 0, iy = 0, iz = 0;
      get_cell_x_y_z(cell_id, ix, iy, iz);

      if ( count != 0 && count % 1000000 == 0)
      {
        krinolog << "Creating local cell " << count << " ." << stk::diag::dendl;
      }
      ++count;

      for (int cell_node_index=0; cell_node_index<num_nodes_per_cell; ++cell_node_index)
      {
        const size_t ax =  ix+hex_cell_node_locations[cell_node_index][0];
        const size_t ay =  iy+hex_cell_node_locations[cell_node_index][1];
        const size_t az =  iz+hex_cell_node_locations[cell_node_index][2];

        cell_node_ids[cell_node_index] = get_node_id(ax, ay, az);
      }
      for (auto && BCC_cell_neighbor : BCC_cell_neighbors)
      {
        cell_node_ids[BCC_cell_neighbor.first] = get_BCC_node_id(ix, iy, iz, BCC_cell_neighbor.second[0], BCC_cell_neighbor.second[1], BCC_cell_neighbor.second[2]);
      }

      if (ix == 0) build_face_tets(cell_id, ix, iy, iz, 0, cell_node_ids);
      if (iy == 0) build_face_tets(cell_id, ix, iy, iz, 2, cell_node_ids);
      if (iz == 0) build_face_tets(cell_id, ix, iy, iz, 4, cell_node_ids);
      build_face_tets(cell_id, ix, iy, iz, 1, cell_node_ids);
      build_face_tets(cell_id, ix, iy, iz, 3, cell_node_ids);
      build_face_tets(cell_id, ix, iy, iz, 5, cell_node_ids);

      for (int cell_node_index=0; cell_node_index<num_nodes_per_cell; ++cell_node_index)
      {
        const size_t ax = ix+hex_cell_node_locations[cell_node_index][0];
        const size_t ay = iy+hex_cell_node_locations[cell_node_index][1];
        const size_t az = iz+hex_cell_node_locations[cell_node_index][2];

        setup_cell_node(ax, ay, az);
      }
      for (auto && BCC_cell_neighbor : BCC_cell_neighbors)
      {
        setup_BCC_node(ix, iy, iz, BCC_cell_neighbor.second[0], BCC_cell_neighbor.second[1], BCC_cell_neighbor.second[2]);
      }
    }
    stk::tools::fix_node_sharing_via_search(*m_mesh);
  }
  m_mesh->modification_end();

  std::vector<size_t> counts;
  stk::mesh::count_entities( m_meta->locally_owned_part(), *m_mesh, counts );
  krinolog << "Generated mesh with " << counts[stk::topology::ELEM_RANK] << " local elements and " << counts[stk::topology::NODE_RANK] << " local nodes." << stk::diag::dendl;
}

void
BoundingBoxMesh::populate_cell_based_mesh()
{ /* %TRACE[ON]% */ Trace trace__("krino::BoundingBoxMesh::populate_cell_based_mesh()"); /* %TRACE% */
  STK_ThrowRequire(m_mesh);
  set_is_cell_edge_function_for_cell_based_mesh();

  const int p_size = m_mesh->parallel_size();
  const int p_rank = m_mesh->parallel_rank();
  const int dim = m_element_topology.dimension();
  stk::mesh::FieldBase const* coord_field = m_meta->coordinate_field();

  std::vector<std::array<int,3>> hex_cell_node_locations = { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1} };
  std::vector<std::vector<int>> hex_cell_elem_nodes = {{0, 1, 2, 3, 4, 5, 6, 7}};

  std::vector<std::vector<int>> tet_even_cell_elem_nodes = {{0, 1, 2, 5},
        {0, 2, 7, 5},
        {0, 2, 3, 7},
        {0, 5, 7, 4},
        {2, 7, 5, 6}};
  std::vector<std::vector<int>> tet_odd_cell_elem_nodes = {{0, 1, 3, 4},
        {1, 2, 3, 6},
        {1, 3, 4, 6},
        {3, 4, 6, 7},
        {1, 6, 4, 5}};

  std::vector<std::array<int,3>> quad_cell_node_locations = { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0} };
  std::vector<std::vector<int>> quad_cell_elem_nodes = {{0, 1, 2, 3}};

  std::vector<std::vector<int>> tri_even_cell_elem_nodes = {{0, 1, 2}, {0, 2, 3}};
  std::vector<std::vector<int>> tri_odd_cell_elem_nodes = {{0, 1, 3}, {1, 2, 3}};

  const std::vector<std::array<int,3>> & cell_node_locations = (dim == 2) ? quad_cell_node_locations : hex_cell_node_locations;
  const std::vector<std::vector<int>> & even_cell_elem_nodes =
      (m_element_topology == stk::topology::TRIANGLE_3_2D) ? tri_even_cell_elem_nodes :
      ((m_element_topology == stk::topology::QUADRILATERAL_4_2D) ? quad_cell_elem_nodes :
      ((m_element_topology == stk::topology::TETRAHEDRON_4) ? tet_even_cell_elem_nodes :
      hex_cell_elem_nodes));
  const std::vector<std::vector<int>> & odd_cell_elem_nodes =
      (m_element_topology == stk::topology::TRIANGLE_3_2D) ? tri_odd_cell_elem_nodes :
      ((m_element_topology == stk::topology::QUADRILATERAL_4_2D) ? quad_cell_elem_nodes :
      ((m_element_topology == stk::topology::TETRAHEDRON_4) ? tet_odd_cell_elem_nodes :
      hex_cell_elem_nodes));

  std::unordered_map<stk::mesh::EntityId, std::vector<int>> nodes_to_procs;
  generate_node_to_processor_map(p_size, p_rank, cell_node_locations, nodes_to_procs);
  const int num_elem_per_cell = even_cell_elem_nodes.size();
  const int num_nodes_per_cell = cell_node_locations.size();
  const int num_nodes_per_elem = even_cell_elem_nodes[0].size();
  const stk::mesh::EntityId elem_id_start = 1;
  std::vector<stk::mesh::EntityId> cell_node_ids(num_nodes_per_cell);
  std::vector<stk::mesh::EntityId> elem_nodes(num_nodes_per_elem);

  my_coord_mapping = std::make_unique<CartesianCoordinateMapping>(m_nx, m_ny, m_nz, m_mesh_bbox);

  std::pair<size_t,size_t> proc_cell_range = determine_processor_cells(p_size, p_rank);
  const size_t num_local_cells = proc_cell_range.second-proc_cell_range.first;
  krinolog << "BoundingBoxMesh creating " << num_local_cells << " local cells." << stk::diag::dendl;

  m_mesh->modification_begin();

  {
    size_t count = 0;
    for (size_t cell_id=proc_cell_range.first; cell_id!=proc_cell_range.second; ++cell_id)
    {
      size_t ix = 0, iy = 0, iz = 0;
      get_cell_x_y_z(cell_id, ix, iy, iz);

      const std::vector<std::vector<int>> & cell_elem_nodes = ((ix+iy+iz)%2 == 0) ? even_cell_elem_nodes : odd_cell_elem_nodes;

      if ( count % 1000000 == 0)
      {
        krinolog << "Creating local cell " << count << " ." << stk::diag::dendl;
      }
      ++count;

      for (int cell_node_index=0; cell_node_index<num_nodes_per_cell; ++cell_node_index)
      {
        const size_t ax =  ix+cell_node_locations[cell_node_index][0];
        const size_t ay =  iy+cell_node_locations[cell_node_index][1];
        const size_t az =  (dim == 3) ?  iz + cell_node_locations[cell_node_index][2] : iz;

        cell_node_ids[cell_node_index] = get_node_id(ax, ay, az);
      }

      for (int cell_elem_index=0; cell_elem_index<num_elem_per_cell; ++cell_elem_index)
      {
        stk::mesh::EntityId elem_id = num_elem_per_cell*cell_id + cell_elem_index + elem_id_start;
        for (int elem_node_index = 0; elem_node_index<num_nodes_per_elem; ++elem_node_index)
        {
          elem_nodes[elem_node_index] = cell_node_ids[cell_elem_nodes[cell_elem_index][elem_node_index]];
        }
        stk::mesh::declare_element( *m_mesh, m_elem_parts, elem_id, elem_nodes );

        for (auto && node_id : elem_nodes)
        {
          stk::mesh::Entity const node = m_mesh->get_entity( stk::topology::NODE_RANK, node_id );
          m_mesh->change_entity_parts(node, m_node_parts);

          auto map_it = nodes_to_procs.find(node_id);
          if (map_it != nodes_to_procs.end())
          {
            for (auto other_proc : map_it->second)
            {
              m_mesh->add_node_sharing(node, other_proc);
            }
          }

          // Compute and assign coordinates to the node
          get_node_x_y_z(node_id, ix, iy, iz);

          double * coord_data = static_cast<double*>(stk::mesh::field_data( *coord_field, node ));
          my_coord_mapping->get_node_coordinates(coord_data, dim, {ix, iy, iz});
        }
      }
    }
  }
  m_mesh->modification_end();

  std::vector<size_t> counts;
  stk::mesh::count_entities( m_meta->locally_owned_part(), *m_mesh, counts );
  krinolog << "Generated mesh with " << counts[stk::topology::ELEM_RANK] << " local elements and " << counts[stk::topology::NODE_RANK] << " local nodes." << stk::diag::dendl;
}

void
BoundingBoxMesh::declare_domain_side_parts(const stk::mesh::Part & blockPart)
{
  AuxMetaData & aux_meta = AuxMetaData::get(*m_meta);

  stk::topology side_topology = m_element_topology.side_topology();

  mySideParts.clear();
  mySideParts.push_back(&aux_meta.declare_io_part_with_topology("Xminus", side_topology));
  mySideParts.push_back(&aux_meta.declare_io_part_with_topology("Xplus", side_topology));
  mySideParts.push_back(&aux_meta.declare_io_part_with_topology("Yminus", side_topology));
  mySideParts.push_back(&aux_meta.declare_io_part_with_topology("Yplus", side_topology));
  if (m_meta->spatial_dimension() == 3)
  {
    mySideParts.push_back(&aux_meta.declare_io_part_with_topology("Zminus", side_topology));
    mySideParts.push_back(&aux_meta.declare_io_part_with_topology("Zplus", side_topology));
  }

  for (auto && sidePart : mySideParts)
    m_meta->set_surface_to_block_mapping(sidePart, {&blockPart});
}

void BoundingBoxMesh::require_has_flat_boundaries() const
{
  STK_ThrowRequireMsg(has_flat_boundaries(), "Domain sides can only be added for CUBIC or FLAT_WALLED_BC generated meshes.");
}

static bool equal_within_tol(const double x1, const double x2, const double tol)
{
  return std::abs(x1-x2) < tol;
}

void
BoundingBoxMesh::create_domain_sides()
{
  if (mySideParts.empty())
    return;

  require_has_flat_boundaries();

  stk::topology side_topology = m_element_topology.side_topology();
  STK_ThrowRequire(mySideParts.size() >= m_meta->spatial_dimension()*2);

  std::vector<stk::mesh::Entity> sides;
  stk::mesh::get_selected_entities( AuxMetaData::get(*m_meta).exposed_boundary_part() & m_meta->locally_owned_part(), m_mesh->buckets( m_meta->side_rank() ), sides );
  std::vector<stk::mesh::PartVector> add_parts(sides.size());
  std::vector<stk::mesh::PartVector> remove_parts(sides.size());
  stk::mesh::FieldBase const* coord_field = m_meta->coordinate_field();

  const auto & min = m_mesh_bbox.get_min();
  const auto & max = m_mesh_bbox.get_max();
  const double relativeTol = 0.01;  // should not be at all sensitive to this tolerance
  const std::array<double,3> tol {(max[0]-min[0])/m_nx*relativeTol, (max[1]-min[1])/m_ny*relativeTol, ((m_meta->spatial_dimension() == 3) ? ((max[2]-min[2])/m_nz*relativeTol) : 0.)};

  const unsigned num_side_nodes = side_topology.num_nodes();
  for (size_t iside=0; iside<sides.size(); ++iside)
  {
    // determine which side of domain this side is on
    const stk::mesh::Entity * side_nodes = m_mesh->begin_nodes(sides[iside]);
    std::array<int,3> domainSide {-2,-2,-2};
    for (unsigned n=0; n<num_side_nodes; ++n)
    {
      double * coords = static_cast<double*>(stk::mesh::field_data( *coord_field, side_nodes[n] ));
      for (unsigned i=0; i<m_meta->spatial_dimension(); ++i )
      {
        if (equal_within_tol(coords[i],min[i],tol[i]))
        {
          STK_ThrowRequire(domainSide[i] != 1);
          if (domainSide[i] == -2) domainSide[i] = -1;
        }
        else if (equal_within_tol(coords[i],max[i],tol[i]))
        {
          STK_ThrowRequire(domainSide[i] != -1);
          if (domainSide[i] == -2) domainSide[i] = 1;
        }
        else
        {
          domainSide[i] = 0;
        }
      }
    }

    const int num_sides_set = ((domainSide[0]==-1||domainSide[0]==1) ? 1 : 0) + ((domainSide[1]==-1||domainSide[1]==1) ? 1 : 0) + ((domainSide[2]==-1||domainSide[2]==1) ? 1 : 0);
    STK_ThrowRequire(num_sides_set == 1);
    stk::mesh::Part * side_part = nullptr;
    if (domainSide[0] == -1) side_part = mySideParts[0];
    else if (domainSide[0] ==  1) side_part = mySideParts[1];
    else if (domainSide[1] == -1) side_part = mySideParts[2];
    else if (domainSide[1] ==  1) side_part = mySideParts[3];
    else if (domainSide[2] == -1) side_part = mySideParts[4];
    else if (domainSide[2] ==  1) side_part = mySideParts[5];

    add_parts[iside].push_back(side_part);
  }

  m_mesh->batch_change_entity_parts(sides, add_parts, remove_parts);
}

void
BoundingBoxMesh::get_cell_x_y_z( stk::mesh::EntityId cell_id, size_t &ix , size_t &iy , size_t &iz ) const
{
  ix = cell_id % m_nx;
  cell_id /= m_nx;

  iy = cell_id % m_ny;
  cell_id /= m_ny;

  iz = cell_id;
}

std::array<size_t,3>
BoundingBoxMesh::get_node_x_y_z( stk::mesh::EntityId entity_id, const std::array<size_t,3> & N )
{
  const size_t node_id_start = 1;
  entity_id -= node_id_start;

  std::array<size_t,3> indices;
  indices[0] = entity_id % (N[0]+1);
  entity_id /= (N[0]+1);

  indices[1] = entity_id % (N[1]+1);
  entity_id /= (N[2]+1);

  indices[2] = entity_id;
  return indices;
}

void
BoundingBoxMesh::get_node_x_y_z( stk::mesh::EntityId entity_id, size_t &ix , size_t &iy , size_t &iz ) const
{
  auto indices = get_node_x_y_z(entity_id, {{m_nx,m_ny,m_nz}});
  ix = indices[0];
  iy = indices[1];
  iz = indices[2];
}

std::pair<size_t, size_t>
BoundingBoxMesh::determine_processor_cells(const int p_size, const int p_rank) const
{
  const size_t Ntot = m_nx*m_ny*m_nz;
  const size_t beg_elem = ( Ntot * p_rank ) / p_size ;
  const size_t end_elem = ( Ntot * ( p_rank + 1 ) ) / p_size ;
  return std::make_pair(beg_elem, end_elem);
}

void
BoundingBoxMesh::generate_node_to_processor_map(const int p_size,
    const int p_rank,
    const std::vector<std::array<int,3>> & cell_node_locations,
    std::unordered_map<stk::mesh::EntityId,
    std::vector<int>> & nodes_to_procs) const
{
  std::unordered_set<stk::mesh::EntityId> locally_used_nodes;

  std::pair<size_t,size_t> proc_cell_range = determine_processor_cells(p_size, p_rank);

  // First create set of all nodes used by local elements
  for (size_t cell_id=proc_cell_range.first; cell_id!=proc_cell_range.second; ++cell_id)
  {
    size_t ix = 0, iy = 0, iz = 0;
    get_cell_x_y_z(cell_id, ix, iy, iz);

    for (size_t cell_node_index=0; cell_node_index<cell_node_locations.size(); ++cell_node_index)
    {
      stk::mesh::EntityId nodeId = get_node_id( ix+cell_node_locations[cell_node_index][0], iy+cell_node_locations[cell_node_index][1] , iz+cell_node_locations[cell_node_index][2] );
      locally_used_nodes.insert(nodeId);
    }
  }

  // Now add entries for nodes that are also used on other procs
  for (int other_p_rank = 0; other_p_rank < p_size; ++other_p_rank)
  {
    if (other_p_rank == p_rank) continue;

    std::pair<size_t,size_t> other_proc_cell_range = determine_processor_cells(p_size, other_p_rank);

    for (size_t cell_id=other_proc_cell_range.first; cell_id!=other_proc_cell_range.second; ++cell_id)
    {
      size_t ix = 0, iy = 0, iz = 0;
      get_cell_x_y_z(cell_id, ix, iy, iz);

      for (size_t cell_node_index=0; cell_node_index<cell_node_locations.size(); ++cell_node_index)
      {
        stk::mesh::EntityId nodeId = get_node_id( ix+cell_node_locations[cell_node_index][0], iy+cell_node_locations[cell_node_index][1] , iz+cell_node_locations[cell_node_index][2] );
        if (locally_used_nodes.find(nodeId) != locally_used_nodes.end())
        {
          std::vector<int> & node_procs = nodes_to_procs[nodeId];
          if (std::find(node_procs.begin(), node_procs.end(), other_p_rank) == node_procs.end())
          {
            node_procs.push_back(other_p_rank);
          }
        }
      }
    }
  }
}

} // namespace krino
