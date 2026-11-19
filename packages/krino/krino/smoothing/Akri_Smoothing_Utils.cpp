// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Smoothing_Utils.hpp>

#include <array>
#include <tuple>
#include <vector>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_MeshHelpers.hpp>
#include <stk_mesh/base/EntityLess.hpp>

namespace krino {

void fill_elem_node_locations(const unsigned spatialDim, const DistributedVector& nodeLocs, const ElemNodeIndices & elemNodes, std::vector<const double*> & elemNodeLocs)
{
  elemNodeLocs.resize(elemNodes.size());
  unsigned i=0;
  for (auto elemNode : elemNodes)
    elemNodeLocs[i++] = &nodeLocs[size_t(spatialDim) * elemNode];
}

void fill_elem_node_locations(const unsigned spatialDim, const DistributedVector& nodeLocs, const MeshElemNodeIndices & meshElemNodeIndices, const unsigned elemIndex, std::vector<const double*> & elemNodeLocs)
{
  const ElemNodeIndices elemNodeIndices = meshElemNodeIndices.element_node_indices(elemIndex);
  fill_elem_node_locations(spatialDim, nodeLocs, elemNodeIndices, elemNodeLocs);
}

static unsigned get_index_in_sorted_entities(const stk::mesh::BulkData &mesh,
  const std::vector<stk::mesh::Entity>::const_iterator sortedEntitiesBeg,
  const std::vector<stk::mesh::Entity>::const_iterator sortedEntitiesEnd,
  const stk::mesh::Entity entity)
{
  auto iter = std::lower_bound(sortedEntitiesBeg, sortedEntitiesEnd, entity, stk::mesh::EntityLess(mesh));
  STK_ThrowAssert(iter != sortedEntitiesEnd && *iter == entity);
  return std::distance(sortedEntitiesBeg, iter);
}

unsigned get_index_of_node_in_sorted_owned_or_shared_nodes(const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    const stk::mesh::Entity node)
{
  return (mesh.bucket(node).owned()) ?
      get_index_in_sorted_entities(mesh, nodes.begin_local(), nodes.end_local(), node) :
      (nodes.local_size() + get_index_in_sorted_entities(mesh, nodes.begin_remote(), nodes.end_remote(), node));
}

DistributedEntityVector get_sorted_owned_nodes_and_unowned_shared_nodes(const stk::mesh::BulkData &mesh, const stk::mesh::Selector & elemSelector)
{
  const bool doSortById = true;
  const stk::mesh::Selector ownedSelector = elemSelector & mesh.mesh_meta_data().locally_owned_part();
  std::vector<stk::mesh::Entity> ownedNodes;
  stk::mesh::get_entities(mesh, stk::topology::NODE_RANK, ownedSelector, ownedNodes, doSortById);
  const stk::mesh::Selector sharedUnownedSelector = elemSelector & mesh.mesh_meta_data().globally_shared_part() & !mesh.mesh_meta_data().locally_owned_part();
  std::vector<stk::mesh::Entity> sharedUnownedNodes;
  stk::mesh::get_entities(mesh, stk::topology::NODE_RANK, sharedUnownedSelector, sharedUnownedNodes, doSortById);

  DistributedEntityVector nodes(mesh.parallel(), ownedNodes, sharedUnownedNodes);
  return nodes;
}

bool is_vertex_node_of_its_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Entity node)
{
  // In a conforming mesh a node is a vertex of all its elements or of none.
  const StkMeshEntities nodeElems{mesh.begin_elements(node), mesh.end_elements(node)};
  return !nodeElems.empty() &&
      is_vertex_node_of_element(
          mesh.bucket(nodeElems[0]).topology(), mesh.begin_element_ordinals(node)[0]);
}

static std::vector<stk::mesh::Entity>
get_vertex_nodes(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Entity> & nodes)
{
  std::vector<stk::mesh::Entity> vertexNodes;
  vertexNodes.reserve(nodes.size());
  for (auto node : nodes)
    if (is_vertex_node_of_its_elements(mesh, node)) vertexNodes.push_back(node);
  return vertexNodes;
}

DistributedEntityVector get_sorted_owned_and_unowned_shared_vertex_nodes(
    const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elemSelector)
{
  const DistributedEntityVector nodes =
      get_sorted_owned_nodes_and_unowned_shared_nodes(mesh, elemSelector);
  const std::vector<stk::mesh::Entity> ownedNodes(nodes.begin_local(), nodes.end_local());
  const std::vector<stk::mesh::Entity> sharedUnownedNodes(nodes.begin_remote(), nodes.end_remote());
  return DistributedEntityVector(mesh.parallel(),
      get_vertex_nodes(mesh, ownedNodes),
      get_vertex_nodes(mesh, sharedUnownedNodes));
}

StkMeshEntities element_vertex_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::Entity elem)
{
  const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);
  return StkMeshEntities{elemNodes, elemNodes + mesh.bucket(elem).topology().base().num_nodes()};
}

NodeLocations get_vertex_node_locations_of_quadratic_elements(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const std::vector<stk::mesh::Entity> & elems)
{
  const unsigned dim = coordsField.dim();
  const auto coordsData = coordsField.data<double>();
  NodeLocations vertexNodeLocations;
  for (auto elem : elems)
  {
    const stk::topology topology = mesh.bucket(elem).topology();
    if (topology == topology.base()) continue;
    for (auto vertexNode : element_vertex_nodes(mesh, elem))
    {
      const auto x = coordsData.entity_values(vertexNode);
      stk::math::Vector3d location = stk::math::Vector3d::ZERO;
      for (unsigned j = 0; j < dim; ++j)
        location[j] = x(stk::mesh::ComponentIdx(j));
      vertexNodeLocations.emplace(vertexNode, location);
    }
  }
  return vertexNodeLocations;
}

template <typename COORDS_DATA>
static bool is_vertex_node_moved(const COORDS_DATA & coordsData,
    const unsigned dim,
    const stk::mesh::Entity vertexNode,
    const NodeLocations & vertexNodeLocationsBeforeSmoothing)
{
  const auto iter = vertexNodeLocationsBeforeSmoothing.find(vertexNode);
  STK_ThrowRequireMsg(iter != vertexNodeLocationsBeforeSmoothing.end(),
      "Missing location before smoothing for a vertex node of a quadratic element.");
  const auto x = coordsData.entity_values(vertexNode);
  for (unsigned j = 0; j < dim; ++j)
    if (x(stk::mesh::ComponentIdx(j)) != iter->second[j])
      return true;
  return false;
}

void set_mid_edge_nodes_of_moved_edges_to_edge_midpoints(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const std::vector<stk::mesh::Entity> & elems,
    const NodeLocations & vertexNodeLocationsBeforeSmoothing)
{
  const unsigned dim = coordsField.dim();
  {
    auto coordsData = coordsField.data<double, stk::mesh::ReadWrite>();
    std::array<unsigned, 3> edgeNodeOrdinals;
    for (auto elem : elems)
    {
      const stk::topology topology = mesh.bucket(elem).topology();
      if (topology == topology.base()) continue;
      const stk::mesh::Entity * elemNodes = mesh.begin_nodes(elem);
      for (unsigned e = 0; e < topology.num_edges(); ++e)
      {
        STK_ThrowRequireMsg(topology.edge_topology(e).num_nodes() == 3,
            "Smoothing quadratic elements needs a single mid-edge node on each edge, not so for "
                << topology << ".");
        topology.edge_node_ordinals(e, edgeNodeOrdinals.data());
        const stk::mesh::Entity midNode = elemNodes[edgeNodeOrdinals[2]];
        if (!mesh.bucket(midNode).owned()) continue;
        const stk::mesh::Entity node0 = elemNodes[edgeNodeOrdinals[0]];
        const stk::mesh::Entity node1 = elemNodes[edgeNodeOrdinals[1]];
        if (!is_vertex_node_moved(coordsData, dim, node0, vertexNodeLocationsBeforeSmoothing) &&
            !is_vertex_node_moved(coordsData, dim, node1, vertexNodeLocationsBeforeSmoothing))
          continue;
        const auto x0 = coordsData.entity_values(node0);
        const auto x1 = coordsData.entity_values(node1);
        auto xMid = coordsData.entity_values(midNode);
        for (unsigned j = 0; j < dim; ++j)
          xMid(stk::mesh::ComponentIdx(j)) =
              0.5 * (x0(stk::mesh::ComponentIdx(j)) + x1(stk::mesh::ComponentIdx(j)));
      }
    }
  }
  parallel_sync_fields(mesh, {&coordsField.field()});
}

static void fill_node_indices_for_element(const stk::mesh::BulkData &mesh,
    const DistributedEntityVector & nodes,
    const stk::mesh::Entity elem,
    unsigned * elemNodeIndices)
{
  const StkMeshEntities elemNodes = element_vertex_nodes(mesh, elem);
  for (unsigned i=0; i<elemNodes.size(); ++i)
    elemNodeIndices[i] = get_index_of_node_in_sorted_owned_or_shared_nodes(mesh, nodes, elemNodes[i]);
}

void MeshElemNodeIndices::build_for_mesh(const stk::mesh::BulkData & mesh,
    const std::vector<stk::mesh::Entity> & elems,
    const DistributedEntityVector & nodes)
{
  myElemNodeStartingIndices.reserve(elems.size()+1);
  unsigned numIndices = 0;
  for (auto & elem : elems)
  {
    myElemNodeStartingIndices.push_back(numIndices);
    numIndices += mesh.bucket(elem).topology().base().num_nodes(); // vertex nodes only
  }
  myElemNodeStartingIndices.push_back(numIndices);

  myElemNodeIndices.resize(numIndices);
  for (unsigned i=0; i<elems.size(); ++i)
    fill_node_indices_for_element(mesh, nodes, elems[i], &myElemNodeIndices[myElemNodeStartingIndices[i]]);
}

DistributedVector gather_node_coordinates(const stk::mesh::BulkData &mesh,
    const CoordinatesFieldRef coordsField,
    const DistributedEntityVector & nodes)
{
  const unsigned dim = coordsField.dim();
  DistributedVector nodesCoords(mesh.parallel(), dim*nodes.size(), dim*nodes.local_size());
  auto coordsData = coordsField.data<double>();
  for (unsigned i=0; i<nodes.size(); ++i)
  {
    const stk::mesh::Entity node = nodes[i];
    auto nodeCoordData = coordsData.entity_values(node);
    for (unsigned j=0; j<dim; ++j)
      nodesCoords[dim*i+j] = nodeCoordData(stk::mesh::ComponentIdx(j));
  }
  return nodesCoords;
}

void set_node_coordinates(const stk::mesh::BulkData & mesh,
    const CoordinatesFieldRef coordsField,
    const DistributedEntityVector & nodes,
    const DistributedVector & nodesCoords)
{
  const unsigned dim = coordsField.dim();
  {
    auto coordsData = coordsField.data<double, stk::mesh::ReadWrite>();
    for (unsigned i = 0; i < nodes.local_size(); ++i)
    {
      auto nodeCoordData = coordsData.entity_values(nodes[i]);
      for (unsigned j = 0; j < dim; ++j)
        nodeCoordData(stk::mesh::ComponentIdx(j)) = nodesCoords[dim * i + j];
    }
  }
  parallel_sync_fields(mesh, {&coordsField.field()});
}

static void pack_shared_node_sensitivity_contributions_for_owners(const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    const DistributedVector & sens,
    stk::CommSparse &commSparse)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (size_t i=nodes.local_size(); i<nodes.size(); ++i)
    {
      const double * nodeSens = sens.data() + dim*i;
      const int owner = mesh.parallel_owner_rank(nodes[i]);

      commSparse.send_buffer(owner).pack(mesh.identifier(nodes[i]));
      commSparse.send_buffer(owner).pack(nodeSens, dim);
    }
  });
}

static void unpack_shared_node_sensitivity_contributions(const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    DistributedVector & sens,
    stk::CommSparse &commSparse)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId nodeId;
      commSparse.recv_buffer(procId).unpack(nodeId);
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeId);
      STK_ThrowRequire(mesh.is_valid(node));
      const unsigned index = get_index_in_sorted_entities(mesh, nodes.begin_local(), nodes.end_local(), node);
      std::array<double, 3> contrib;
      commSparse.recv_buffer(procId).unpack(contrib.data(), dim);
      for (unsigned i=0; i<dim; ++i)
        sens[dim*index + i] += contrib[i];
    }
  });
}

static void communicate_shared_node_sensitivity_contributions_to_owners(const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    DistributedVector & sens)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_shared_node_sensitivity_contributions_for_owners(mesh, nodes, sens, commSparse);
  unpack_shared_node_sensitivity_contributions(mesh, nodes, sens, commSparse);
}

static void pack_owned_node_sensitivities_for_sharers(const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    const DistributedVector & sens,
    stk::CommSparse &commSparse)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  std::vector<int> nodeSharedProcs;
  stk::pack_and_communicate(commSparse,[&]()
  {
    for (unsigned i=0; i<nodes.local_size(); ++i)
    {
      if (mesh.bucket(nodes[i]).shared())
      {
        const double * nodeSens = sens.data() + dim*i;
        const stk::mesh::EntityId nodeId = mesh.identifier(nodes[i]);
        mesh.comm_shared_procs(nodes[i], nodeSharedProcs);
        for (int procId : nodeSharedProcs)
        {
          commSparse.send_buffer(procId).pack(nodeId);
          commSparse.send_buffer(procId).pack(nodeSens, dim);
        }
      }
    }
  });
}

static void unpack_shared_node_sensitivities(const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    DistributedVector & sens,
    stk::CommSparse &commSparse)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  stk::unpack_communications(commSparse, [&](int procId)
  {
    stk::CommBuffer & buffer = commSparse.recv_buffer(procId);

    while ( buffer.remaining() )
    {
      stk::mesh::EntityId nodeId;
      commSparse.recv_buffer(procId).unpack(nodeId);
      stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, nodeId);
      STK_ThrowRequire(mesh.is_valid(node));
      const unsigned index = get_index_in_sorted_entities(mesh, nodes.begin_remote(), nodes.end_remote(), node);
      double * nodeSens = sens.data() + dim*(nodes.local_size()+index);
      commSparse.recv_buffer(procId).unpack(nodeSens, dim);
    }
  });
}

static void communicate_owned_node_sensitivities_to_sharers(const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    DistributedVector & sens)
{
  stk::CommSparse commSparse(mesh.parallel());
  pack_owned_node_sensitivities_for_sharers(mesh, nodes, sens, commSparse);
  unpack_shared_node_sensitivities(mesh, nodes, sens, commSparse);
}

void communicate_gradient(const stk::mesh::BulkData & mesh,
    const DistributedEntityVector & nodes,
    DistributedVector & sens)
{
  if (mesh.parallel_size() > 1)
  {
    communicate_shared_node_sensitivity_contributions_to_owners(mesh, nodes, sens);
    communicate_owned_node_sensitivities_to_sharers(mesh, nodes, sens);
  }
}

void apply_gradient_filter(const unsigned dim,
    const GradientFilter & gradientFilter,
    const DistributedEntityVector & nodes,
    const DistributedVector & x,
    DistributedVector & g)
{
  if (gradientFilter)
  {
    for (size_t i = 0; i < nodes.size(); ++i)
    {
      stk::math::Vector3d coords(&x[dim * i], dim);
      stk::math::Vector3d dir(&g[dim * i], dim);
      gradientFilter(nodes[i], coords, dir);
      for (unsigned j = 0; j < dim; ++j)
        g[dim * i + j] = dir[j];
    }
  }
}

void apply_projection(const unsigned dim,
    const SolutionProjector & solutionProjector,
    const DistributedEntityVector & nodes,
    DistributedVector & x)
{
  if (solutionProjector)
  {
    for (size_t i = 0; i < nodes.size(); ++i)
    {
      stk::math::Vector3d coords(&x[dim * i], dim);
      solutionProjector(nodes[i], coords);
      for (unsigned j = 0; j < dim; ++j)
        x[dim * i + j] = coords[j];
    }
  }
}


}


